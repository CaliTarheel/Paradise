// Stepping the history: a stage pole per plate, preview, commit, go back, export.
'use strict';

const KIND_COLOURS = {
  ridge: '#d0362b', trench: '#2a5fae', transform: '#2f8a4f', converging: '#c78a12', collision: '#7a4a1c',
};

const stepper = {
  state: null,
  preview: null,
  selected: null,
  mapImage: new Image(),
  plateOf: null,

  async load() {
    const r = await fetch('/api/state');
    const data = await r.json();
    if (data.empty) return;
    this.show(data);
  },

  show(state) {
    this.state = state;
    this.preview = null;
    this.plateOf = fromBase64(state.plate_of);
    this.mapImage.onload = () => this.draw();
    this.mapImage.src = state.map;
    document.getElementById('clock').textContent =
      `${state.name}: ${fmt(state.age_ma)} Ma (step ${state.steps.length - 1})`;
    this.table();
    this.history();
    this.notes(state.warnings, state.events);
  },

  table() {
    const body = document.querySelector('#poles tbody');
    body.innerHTML = '';
    for (const p of this.state.plates) {
      const tr = document.createElement('tr');
      tr.dataset.plate = p.id;
      if (p.id === this.selected) tr.classList.add('selected');
      tr.innerHTML = `<td title="${p.area_mkm2} million km²">${p.id}</td>` +
        ['lat', 'lon', 'rate'].map((k, i) =>
          `<td><input id="pole-${p.id}-${k}" type="number" step="${k === 'rate' ? 0.05 : 1}" value="${p.pole[i]}" aria-label="Plate ${p.id} pole ${k}"></td>`).join('') +
        `<td class="speed" id="speed-${p.id}"></td>`;
      tr.addEventListener('click', () => this.select(p.id));
      tr.querySelectorAll('input').forEach((input) => input.addEventListener('input', () => this.speeds()));
      body.appendChild(tr);
    }
    this.speeds();
  },

  select(plate) {
    this.selected = plate;
    document.querySelectorAll('#poles tr').forEach((tr) => tr.classList.toggle('selected', Number(tr.dataset.plate) === plate));
    this.draw();
  },

  poles() {
    const poles = {};
    for (const p of this.state.plates) {
      const v = ['lat', 'lon', 'rate'].map((k) => Number(document.getElementById(`pole-${p.id}-${k}`).value));
      poles[p.id] = v;
    }
    return poles;
  },

  // Fastest speed on each plate's outline for its pole, cm/yr (1 degree/Myr is 11.1 cm/yr at 90 degrees).
  speeds() {
    const poles = this.poles();
    for (const p of this.state.plates) {
      const [lat, lon, rate] = poles[p.id];
      const axis = unit(lon, lat);
      let fastest = 0;
      for (const o of this.state.lines.outlines) {
        if (o.plate !== p.id) continue;
        for (const [x, y] of o.points) {
          const q = unit(x, y);
          const s = Math.sqrt(Math.max(0, 1 - dot(axis, q) ** 2));
          fastest = Math.max(fastest, s);
        }
      }
      const cm = Math.abs(rate) * Math.PI / 180 * 6371 * fastest / 10;
      const cell = document.getElementById(`speed-${p.id}`);
      cell.textContent = cm.toFixed(1);
      cell.classList.toggle('fast', cm > 20);
    }
  },

  history() {
    const list = document.getElementById('steps');
    list.innerHTML = '';
    for (const s of this.state.steps) {
      const li = document.createElement('li');
      li.textContent = `${fmt(s.age_ma)} Ma ${s.note ? '— ' + s.note : ''} `;
      if (s.index < this.state.steps.length - 1) {
        const b = document.createElement('button');
        b.textContent = 'Go back here';
        b.addEventListener('click', async () => {
          const result = await api('/api/back', { step: s.index }, 'Going back…');
          if (result) this.show(result);
        });
        li.appendChild(b);
      }
      list.appendChild(li);
    }
  },

  notes(warnings, events) {
    const box = document.getElementById('warnings');
    box.innerHTML = '';
    for (const w of warnings || []) {
      const p = document.createElement('p'); p.textContent = w.text; box.appendChild(p);
    }
    const seen = new Set();
    for (const e of events || []) {
      const text = describe(e);
      if (!text || seen.has(text)) continue;
      seen.add(text);
      const p = document.createElement('p'); p.textContent = text; box.appendChild(p);
    }
  },

  draw() {
    const canvas = document.getElementById('map');
    const ctx = canvas.getContext('2d');
    const w = canvas.width, h = canvas.height;
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(this.mapImage, 0, 0, w, h);
    const lines = this.state.lines;
    for (const o of lines.outlines) {
      const selected = o.plate === this.selected;
      polyline(ctx, o.points, w, h, selected ? '#ffffff' : 'rgba(0,0,0,0.35)', selected ? 2.5 : 1);
    }
    for (const c of lines.cratons || []) polyline(ctx, c.points, w, h, 'rgba(70,45,20,0.9)', 1.5, [3, 3]);
    for (const b of lines.boundaries) polyline(ctx, b.points, w, h, KIND_COLOURS[b.kind] || '#000', 2.5);
    if (this.preview) {
      for (const o of this.preview.lines.outlines) polyline(ctx, o.points, w, h, 'rgba(255,255,255,0.9)', 1.5, [6, 4]);
    }
    // The selected plate's pole.
    if (this.selected != null) {
      const [lat, lon] = this.poles()[this.selected];
      const x = (lon + 180) / 360 * w, y = (90 - lat) / 180 * h;
      ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(x, y, 7, 0, 2 * Math.PI); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x - 11, y); ctx.lineTo(x + 11, y); ctx.moveTo(x, y - 11); ctx.lineTo(x, y + 11); ctx.stroke();
    }
    ctx.fillStyle = '#ffffff'; ctx.font = 'bold 14px system-ui'; ctx.textAlign = 'center';
    ctx.shadowColor = '#000'; ctx.shadowBlur = 3;
    for (const p of this.state.plates) {
      ctx.fillText(String(p.id), (p.centre[0] + 180) / 360 * w, (90 - p.centre[1]) / 180 * h);
    }
    ctx.shadowBlur = 0;
  },
};

function describe(e) {
  const pair = e.plates ? `plates ${e.plates[0]} and ${e.plates[1]}` : '';
  switch (e.kind) {
    case 'subduction_initiated': return `A trench formed between ${pair}; plate ${e.overriding} overrides.`;
    case 'collision': return `Continents collided between ${pair}.`;
    case 'suture': return `Plates ${e.plates[0]} and ${e.plates[1]} welded along a suture into plate ${e.into}.`;
    case 'plate_split': return `Plate ${e.plate} came apart; plate ${e.new} is the new piece.`;
    case 'absorbed': return `A piece of plate ${e.plate} (${Math.round(e.area_km2 / 1000)} thousand km²) joined plate ${e.into}.`;
    case 'plate_gone': return `Plate ${e.plate} is gone.`;
    case 'rift': return `Plate ${e.plate} rifted; plate ${e.new} is the new piece.`;
    default: return e.kind ? `${e.kind.replace(/_/g, ' ')}${pair ? ' (' + pair + ')' : ''}` : '';
  }
}

function polyline(ctx, points, w, h, colour, width, dash) {
  ctx.strokeStyle = colour; ctx.lineWidth = width; ctx.setLineDash(dash || []);
  ctx.beginPath();
  let px = null;
  for (const [lon, lat] of points) {
    const x = (lon + 180) / 360 * w, y = (90 - lat) / 180 * h;
    if (px === null || Math.abs(x - px) > w / 2) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    px = x;
  }
  ctx.stroke();
  ctx.setLineDash([]);
}

function unit(lon, lat) {
  const a = lon * Math.PI / 180, b = lat * Math.PI / 180;
  return [Math.cos(b) * Math.cos(a), Math.cos(b) * Math.sin(a), Math.sin(b)];
}
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function fmt(age) { return Math.round(age * 10) / 10; }

document.getElementById('map').addEventListener('click', (e) => {
  if (!stepper.state) return;
  const canvas = e.target, r = canvas.getBoundingClientRect();
  const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
  if (e.shiftKey && stepper.selected != null) {
    document.getElementById(`pole-${stepper.selected}-lat`).value = Math.round(90 - fy * 180);
    document.getElementById(`pole-${stepper.selected}-lon`).value = Math.round(fx * 360 - 180);
    stepper.speeds(); stepper.draw();
    return;
  }
  const width = stepper.state.map_width, height = width / 2;
  const plate = stepper.plateOf[Math.floor(fy * height) * width + Math.floor(fx * width)];
  if (plate) stepper.select(plate);
});

document.getElementById('preview').addEventListener('click', async () => {
  const result = await api('/api/preview', { poles: stepper.poles() }, 'Previewing the next step…');
  if (result) { stepper.preview = result; stepper.draw(); stepper.notes(result.warnings, result.events); }
});

document.getElementById('commit').addEventListener('click', async () => {
  const note = document.getElementById('note').value;
  const result = await api('/api/commit', { poles: stepper.poles(), note }, 'Stepping 10 Myr…');
  if (result) { document.getElementById('note').value = ''; stepper.show(result); }
});

document.getElementById('export').addEventListener('click', async () => {
  const result = await api('/api/export', {}, 'Writing the GPlates project…');
  if (result) {
    document.getElementById('export-result').textContent =
      `Wrote ${result.files.join(', ')} to ${result.directory}. In GPlates, open all of them together.`;
  }
});

stepper.load().then(() => { if (stepper.state) showView('step'); });
