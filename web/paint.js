// Painting the starting world on an equirectangular grid of 0.5 degree cells.
// Rows run from the north pole down; columns from 180 W eastwards, wrapping round.
'use strict';

const W = 720, H = 360, RADIUS_KM = 6371;
const OCEAN = 0, CONTINENT = 1, CRATON = 2;

const paint = {
  plates: new Uint16Array(W * H).fill(1),
  crust: new Uint16Array(W * H),
  plate: 1,
  plateCount: 1,
  undo: [],
};

const canvas = document.getElementById('paint');
const ctx = canvas.getContext('2d');
const image = ctx.createImageData(W, H);

function plateColour(p) {
  // Golden-angle hues, so neighbouring numbers differ.
  const hue = (p * 137.508) % 360;
  return hslToRgb(hue / 360, 0.45, 0.55);
}

function hslToRgb(h, s, l) {
  const f = (n) => {
    const k = (n + h * 12) % 12;
    const a = s * Math.min(l, 1 - l);
    return Math.round(255 * (l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1))));
  };
  return [f(0), f(8), f(4)];
}

function render() {
  const d = image.data;
  for (let i = 0; i < W * H; i++) {
    const [r, g, b] = plateColour(paint.plates[i]);
    const c = paint.crust[i];
    let k;
    if (c === OCEAN) k = [r * 0.45 + 40, g * 0.45 + 70, b * 0.45 + 110];
    else if (c === CONTINENT) k = [r * 0.55 + 105, g * 0.55 + 95, b * 0.55 + 60];
    else k = [r * 0.35 + 60, g * 0.35 + 50, b * 0.35 + 40];
    d[i * 4] = k[0]; d[i * 4 + 1] = k[1]; d[i * 4 + 2] = k[2]; d[i * 4 + 3] = 255;
  }
  // Outline plate edges so the plates read through the crust colours.
  for (let y = 0; y < H; y++) {
    for (let x = 0; x < W; x++) {
      const i = y * W + x, p = paint.plates[i];
      const right = paint.plates[y * W + (x + 1) % W];
      const down = y + 1 < H ? paint.plates[i + W] : p;
      if (p !== right || p !== down) { d[i * 4] = d[i * 4 + 1] = d[i * 4 + 2] = 20; }
    }
  }
  ctx.putImageData(image, 0, 0);
}

function cellCentre(x, y) {
  return [(-180 + (x + 0.5) * 360 / W) * Math.PI / 180, (90 - (y + 0.5) * 180 / H) * Math.PI / 180];
}

function arcKm(lon1, lat1, lon2, lat2) {
  const s = Math.sin((lat2 - lat1) / 2) ** 2 + Math.cos(lat1) * Math.cos(lat2) * Math.sin((lon2 - lon1) / 2) ** 2;
  return 2 * RADIUS_KM * Math.asin(Math.min(1, Math.sqrt(s)));
}

function layer() { return document.querySelector('input[name=layer]:checked').value; }
function tool() { return document.querySelector('input[name=tool]:checked').value; }
function value() {
  return layer() === 'plates' ? paint.plate : Number(document.querySelector('input[name=crust]:checked').value);
}
function grid() { return layer() === 'plates' ? paint.plates : paint.crust; }

function dab(x, y) {
  const km = Number(document.getElementById('brush-km').value);
  const [lon0, lat0] = cellCentre(x, y);
  const rows = Math.ceil(km / (Math.PI * RADIUS_KM / H)) + 1;
  const g = grid(), v = value();
  for (let yy = Math.max(0, y - rows); yy <= Math.min(H - 1, y + rows); yy++) {
    for (let xx = 0; xx < W; xx++) {
      const [lon, lat] = cellCentre(xx, yy);
      if (arcKm(lon0, lat0, lon, lat) <= km) g[yy * W + xx] = v;
    }
  }
}

function fill(x, y) {
  const g = grid(), v = value(), target = g[y * W + x];
  if (target === v) return;
  const stack = [y * W + x];
  while (stack.length) {
    const i = stack.pop();
    if (g[i] !== target) continue;
    g[i] = v;
    const cx = i % W, cy = (i - cx) / W;
    stack.push(cy * W + (cx + 1) % W, cy * W + (cx + W - 1) % W);
    if (cy > 0) stack.push(i - W);
    if (cy < H - 1) stack.push(i + W);
  }
}

function remember() {
  paint.undo.push([paint.plates.slice(), paint.crust.slice()]);
  if (paint.undo.length > 40) paint.undo.shift();
}

function cellAt(event) {
  const r = canvas.getBoundingClientRect();
  return [Math.floor((event.clientX - r.left) / r.width * W), Math.floor((event.clientY - r.top) / r.height * H)];
}

let painting = false;
canvas.addEventListener('pointerdown', (e) => {
  remember();
  const [x, y] = cellAt(e);
  if (tool() === 'fill') { fill(x, y); render(); return; }
  painting = true;
  canvas.setPointerCapture(e.pointerId);
  dab(x, y); render();
});
canvas.addEventListener('pointermove', (e) => {
  if (!painting) return;
  const [x, y] = cellAt(e);
  if (x < 0 || y < 0 || x >= W || y >= H) return;
  dab(x, y); render();
});
canvas.addEventListener('pointerup', () => { painting = false; });

function swatches() {
  const box = document.getElementById('plate-swatches');
  box.innerHTML = '';
  for (let p = 1; p <= paint.plateCount; p++) {
    const b = document.createElement('button');
    b.className = 'swatch' + (p === paint.plate ? ' selected' : '');
    const [r, g, bl] = plateColour(p);
    b.style.background = `rgb(${r},${g},${bl})`;
    b.textContent = p;
    b.title = 'Plate ' + p;
    b.addEventListener('click', () => { paint.plate = p; swatches(); });
    box.appendChild(b);
  }
}

document.getElementById('new-plate').addEventListener('click', () => {
  paint.plateCount += 1; paint.plate = paint.plateCount; swatches();
});
document.querySelectorAll('input[name=layer]').forEach((r) => r.addEventListener('change', () => {
  document.getElementById('plate-tools').hidden = layer() !== 'plates';
  document.getElementById('crust-tools').hidden = layer() !== 'crust';
}));
const brush = document.getElementById('brush-km');
brush.addEventListener('input', () => { document.getElementById('brush-km-out').value = brush.value; });
document.getElementById('undo').addEventListener('click', () => {
  const last = paint.undo.pop();
  if (last) { [paint.plates, paint.crust] = last; render(); }
});

// A seeded random number generator, so the same seed gives the same plates.
function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

document.getElementById('voronoi').addEventListener('click', () => {
  remember();
  const n = Math.max(2, Math.min(60, Number(document.getElementById('voronoi-count').value)));
  const rand = mulberry32(Number(document.getElementById('voronoi-seed').value));
  const seeds = [];
  for (let k = 0; k < n; k++) {
    const z = 2 * rand() - 1, t = 2 * Math.PI * rand(), r = Math.sqrt(1 - z * z);
    seeds.push([r * Math.cos(t), r * Math.sin(t), z]);
  }
  for (let y = 0; y < H; y++) {
    for (let x = 0; x < W; x++) {
      const [lon, lat] = cellCentre(x, y);
      const p = [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
      let best = 0, bestDot = -2;
      for (let k = 0; k < n; k++) {
        const dot = p[0] * seeds[k][0] + p[1] * seeds[k][1] + p[2] * seeds[k][2];
        if (dot > bestDot) { bestDot = dot; best = k; }
      }
      paint.plates[y * W + x] = best + 1;
    }
  }
  paint.plateCount = n; paint.plate = 1;
  swatches(); render();
});

function toBase64(array) {
  const bytes = new Uint8Array(array.buffer);
  let s = '';
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
  return btoa(s);
}

function fromBase64(text) {
  const s = atob(text), bytes = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) bytes[i] = s.charCodeAt(i);
  return new Uint16Array(bytes.buffer);
}

document.getElementById('build').addEventListener('click', async () => {
  // Plates are numbered by what is painted, so unused numbers do not become empty plates.
  const body = {
    name: document.getElementById('world-name').value || 'Paradise',
    start_age_ma: Number(document.getElementById('start-age').value),
    ocean_age_myr: Number(document.getElementById('ocean-age').value),
    auto_trenches: document.getElementById('auto-trenches').checked,
    terrane_max_km2: 1000 * Number(document.getElementById('terrane-max').value),
    coast_detail_km: Number(document.getElementById('coast-detail').value),
    width: W, height: H,
    plates: toBase64(paint.plates), crust: toBase64(paint.crust),
  };
  const result = await api('/api/build', body, 'Building the world…');
  if (result) { showView('step'); stepper.show(result); }
});

// Restore the painting of the project the server has open, if any.
async function loadPainting() {
  const r = await fetch('/api/start');
  const data = r.ok ? await r.json() : { empty: true };
  if (!data.empty && data.width === W && data.height === H) {
    paint.plates = fromBase64(data.plates);
    paint.crust = fromBase64(data.crust);
    paint.plateCount = paint.plates.reduce((a, b) => Math.max(a, b), 1);
    document.getElementById('world-name').value = data.name;
    document.getElementById('start-age').value = data.start_age_ma;
    document.getElementById('ocean-age').value = data.ocean_age_myr;
    document.getElementById('auto-trenches').checked = data.auto_trenches !== false;
    if (data.coast_detail_km != null) document.getElementById('coast-detail').value = data.coast_detail_km;
    if (data.terrane_max_km2 != null) document.getElementById('terrane-max').value = Math.round(data.terrane_max_km2 / 1000);
  }
  swatches(); render();
}

function setStatus(text, error) {
  const s = document.getElementById('status');
  s.textContent = text; s.className = error ? 'error' : '';
}

async function api(path, body, busy) {
  setStatus(busy || 'Working…');
  try {
    const r = await fetch(path, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body || {}) });
    const data = await r.json();
    if (!r.ok) { setStatus(data.error || 'The server refused that.', true); return null; }
    setStatus(data.message || 'Done');
    return data;
  } catch (e) {
    setStatus('Could not reach the Paradise server. Is it running?', true);
    return null;
  }
}

function showView(name) {
  document.querySelectorAll('.tab').forEach((t) => t.classList.toggle('active', t.dataset.view === name));
  document.getElementById('view-paint').hidden = name !== 'paint';
  document.getElementById('view-step').hidden = name !== 'step';
}
document.querySelectorAll('.tab').forEach((t) => t.addEventListener('click', () => showView(t.dataset.view)));

loadPainting();
