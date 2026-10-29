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
  // A picked stretch of boundary: the run it lies on, the vertex first clicked, and the stretch [i0, i1].
  pick: null,
  // Declarations for the next step, sent with Preview and Commit.
  declared: [],

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
    this.pick = null;
    this.declared = [];
    this.table();
    this.history();
    this.notes(state.warnings, state.events);
    this.trenchTable(state.trenches);
    this.declaredList();
    this.pickPanel();
  },

  // `preview`: the trenches a preview would leave. They are not where the boundaries are now, so
  // they cannot be picked.
  trenchTable(trenches, preview) {
    const body = document.querySelector('#trenches tbody');
    document.getElementById('trench-mode').textContent = preview ? '(after the preview)'
      : this.state.auto_trenches ? '(they also form by themselves)' : '(only declared ones)';
    body.innerHTML = '';
    for (const t of trenches || []) {
      const tr = document.createElement('tr');
      const cells = [`${t.plates[0]}–${t.plates[1]}`, t.overriding, t.downgoing, t.length_km.toLocaleString(),
        t.closing_cm_yr.toFixed(1), t.slab_age_myr == null ? '—' : t.slab_age_myr.toFixed(0),
        (t.consumed_km2_per_step / 1e6).toFixed(2)];
      for (const c of cells) { const td = document.createElement('td'); td.textContent = c; tr.appendChild(td); }
      if (t.closing_cm_yr < 0) tr.title = 'The plates pull apart across this trench.';
      if (!preview) tr.addEventListener('click', () => {
        const run = { kind: 'trench', plates: t.plates, overriding: t.overriding, points: t.points };
        this.pick = { run, anchor: 0, i0: 0, i1: t.points.length - 1 };
        this.pickPanel(); this.draw();
      });
      body.appendChild(tr);
    }
    if (!body.children.length) {
      const tr = document.createElement('tr');
      tr.innerHTML = '<td colspan="7" class="hint">No trenches yet.</td>';
      body.appendChild(tr);
    }
  },

  // Boundary picking: the first click picks a whole run, a second click on it keeps the stretch between.
  pickAt(x, y) {
    const canvas = document.getElementById('map');
    const w = canvas.width, h = canvas.height;
    let best = null;
    for (const b of this.state.lines.boundaries) {
      for (let i = 0; i < b.points.length; i++) {
        const [px, py] = toCanvas(b.points[i], w, h);
        const d = Math.hypot(px - x, py - y);
        if (!best || d < best.d) best = { d, run: b, i };
      }
      for (let i = 0; i + 1 < b.points.length; i++) {
        const [ax, ay] = toCanvas(b.points[i], w, h), [bx, by] = toCanvas(b.points[i + 1], w, h);
        if (Math.abs(bx - ax) > w / 2) continue;
        const d = segmentDistance(x, y, ax, ay, bx, by);
        if (d < best.d) best = { d, run: b, i: Math.hypot(ax - x, ay - y) < Math.hypot(bx - x, by - y) ? i : i + 1 };
      }
    }
    if (!best || best.d > 14) { this.pick = null; this.pickPanel(); this.draw(); return; }
    if (this.pick && this.pick.run === best.run && this.pick.i0 === 0 && this.pick.i1 === best.run.points.length - 1
        && best.i !== this.pick.anchor) {
      this.pick.i0 = Math.min(this.pick.anchor, best.i);
      this.pick.i1 = Math.max(this.pick.anchor, best.i);
    } else {
      this.pick = { run: best.run, anchor: best.i, i0: 0, i1: best.run.points.length - 1 };
    }
    this.pickPanel(); this.draw();
  },

  pickedPoints() {
    const { run, i0, i1 } = this.pick;
    const a = Math.max(0, i1 === i0 ? i0 - 1 : i0), b = Math.min(run.points.length - 1, i1 === i0 ? i1 + 1 : i1);
    return run.points.slice(a, b + 1);
  },

  pickPanel() {
    const panel = document.getElementById('boundary-pick');
    if (!this.pick) { panel.hidden = true; return; }
    panel.hidden = false;
    const { run } = this.pick;
    const [a, b] = run.plates;
    const km = lengthKm(this.pickedPoints());
    const whole = this.pick.i0 === 0 && this.pick.i1 === run.points.length - 1;
    const state = run.kind === 'trench' ? `trench, plate ${run.overriding} over plate ${run.overriding === a ? b : a}` : run.kind;
    document.getElementById('boundary-text').textContent =
      `Boundary between plates ${a} and ${b} (${state}): ${Math.round(km).toLocaleString()} km picked. ` +
      (whole ? 'Click this boundary again to keep only the stretch between your two clicks.' : 'Click elsewhere to pick again.');
    const left = document.getElementById('subduct-left'), right = document.getElementById('subduct-right');
    left.textContent = `Plate ${a} goes under plate ${b}`;
    right.textContent = `Plate ${b} goes under plate ${a}`;
    left.onclick = () => this.declare('subduction', b);
    right.onclick = () => this.declare('subduction', a);
    const end = document.getElementById('end-subduction');
    end.hidden = run.kind !== 'trench';
    end.onclick = () => this.declare('end_subduction');
    document.getElementById('pick-clear').onclick = () => { this.pick = null; this.pickPanel(); this.draw(); };
  },

  declare(kind, overriding) {
    const event = { kind, plates: this.pick.run.plates.slice(), points: this.pickedPoints() };
    if (kind === 'subduction') event.overriding = overriding;
    this.declared.push(event);
    this.pick = null;
    this.preview = null;
    this.pickPanel(); this.declaredList(); this.trenchTable(this.state.trenches); this.draw();
  },

  declaredList() {
    const list = document.getElementById('declared');
    list.innerHTML = '';
    if (!this.declared.length) {
      list.innerHTML = '<li class="hint">Nothing declared. Pick a boundary on the map to declare a trench.</li>';
      return;
    }
    this.declared.forEach((e, k) => {
      const li = document.createElement('li');
      li.textContent = describeDeclaration(e) + ` (${Math.round(lengthKm(e.points)).toLocaleString()} km)`;
      const b = document.createElement('button');
      b.textContent = 'Remove';
      b.addEventListener('click', () => {
        this.declared.splice(k, 1); this.preview = null;
        this.declaredList(); this.trenchTable(this.state.trenches); this.draw();
      });
      li.appendChild(b);
      list.appendChild(li);
    });
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
      const declared = (s.events || []).map(describeDeclaration).join('; ');
      li.textContent = `${fmt(s.age_ma)} Ma ${s.note ? '— ' + s.note : ''}${declared ? ' [' + declared + ']' : ''} `;
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
    for (const b of lines.boundaries) {
      polyline(ctx, b.points, w, h, KIND_COLOURS[b.kind] || '#000', 2.5);
      if (b.kind === 'trench' && b.overriding != null) teeth(ctx, b.points, w, h, b.overriding === b.plates[0], KIND_COLOURS.trench);
    }
    if (this.preview) {
      for (const o of this.preview.lines.outlines) polyline(ctx, o.points, w, h, 'rgba(255,255,255,0.9)', 1.5, [6, 4]);
      for (const b of this.preview.lines.boundaries) {
        if (b.kind === 'trench' && b.overriding != null) teeth(ctx, b.points, w, h, b.overriding === b.plates[0], 'rgba(255,255,255,0.9)');
      }
    }
    // Declarations waiting for the next step: the stretch, with teeth on the side of the plate that will override.
    for (const e of this.declared) {
      polyline(ctx, e.points, w, h, 'rgba(255,255,255,0.95)', 6);
      polyline(ctx, e.points, w, h, e.kind === 'subduction' ? KIND_COLOURS.trench : '#555', 3, [8, 5]);
      if (e.kind === 'subduction') teeth(ctx, e.points, w, h, e.overriding === e.plates[0], KIND_COLOURS.trench, 11);
    }
    if (this.pick) {
      polyline(ctx, this.pick.run.points, w, h, 'rgba(255,255,255,0.45)', 7);
      polyline(ctx, this.pickedPoints(), w, h, '#ffffff', 7);
      polyline(ctx, this.pickedPoints(), w, h, '#111111', 2);
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

function describeDeclaration(e) {
  const [a, b] = e.plates;
  if (e.kind === 'subduction') return `Trench: plate ${e.overriding === a ? b : a} goes under plate ${e.overriding}`;
  if (e.kind === 'end_subduction') return `End the trench between plates ${a} and ${b}`;
  return e.kind;
}

function toCanvas([lon, lat], w, h) { return [(lon + 180) / 360 * w, (90 - lat) / 180 * h]; }

function segmentDistance(x, y, ax, ay, bx, by) {
  const dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy;
  const t = l2 ? Math.max(0, Math.min(1, ((x - ax) * dx + (y - ay) * dy) / l2)) : 0;
  return Math.hypot(x - ax - t * dx, y - ay - t * dy);
}

function lengthKm(points) {
  let km = 0;
  for (let i = 0; i + 1 < points.length; i++) {
    const a = unit(...points[i]), b = unit(...points[i + 1]);
    km += 6371 * Math.acos(Math.max(-1, Math.min(1, dot(a, b))));
  }
  return km;
}

// Trench teeth on the side of the overriding plate. Walking the points, the boundary's first plate is
// on the left; seen on the map (north up), left of a heading (dx, dy) in canvas pixels is (dy, -dx).
function teeth(ctx, points, w, h, onLeft, colour, size) {
  size = size || 8;
  const spacing = 3 * size;
  ctx.fillStyle = colour;
  let carry = spacing / 2;
  for (let i = 0; i + 1 < points.length; i++) {
    const [ax, ay] = toCanvas(points[i], w, h), [bx, by] = toCanvas(points[i + 1], w, h);
    if (Math.abs(bx - ax) > w / 2) continue;
    const len = Math.hypot(bx - ax, by - ay);
    if (!len) continue;
    const ux = (bx - ax) / len, uy = (by - ay) / len;
    const nx = onLeft ? uy : -uy, ny = onLeft ? -ux : ux;
    let t = carry;
    for (; t < len; t += spacing) {
      const cx = ax + ux * t, cy = ay + uy * t;
      ctx.beginPath();
      ctx.moveTo(cx - ux * size * 0.6, cy - uy * size * 0.6);
      ctx.lineTo(cx + ux * size * 0.6, cy + uy * size * 0.6);
      ctx.lineTo(cx + nx * size, cy + ny * size);
      ctx.closePath();
      ctx.fill();
    }
    carry = t - len;   // how far into the next segment the next tooth goes
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

function mapMode() { return document.querySelector('input[name=map-mode]:checked').value; }

document.querySelectorAll('input[name=map-mode]').forEach((r) => r.addEventListener('change', () => {
  document.getElementById('map-hint').textContent = mapMode() === 'boundary'
    ? 'Click a boundary to pick it. Click it again to keep only the stretch between your two clicks.'
    : "Click a plate to select its row. Shift-click the map to put the selected plate's pole there.";
}));

document.getElementById('map').addEventListener('click', (e) => {
  if (!stepper.state) return;
  const canvas = e.target, r = canvas.getBoundingClientRect();
  const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
  if (mapMode() === 'boundary' && !e.shiftKey) {
    stepper.pickAt(fx * canvas.width, fy * canvas.height);
    return;
  }
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
  const result = await api('/api/preview', { poles: stepper.poles(), events: stepper.declared }, 'Previewing the next step…');
  if (result) {
    stepper.preview = result; stepper.draw();
    stepper.notes(result.warnings, result.events); stepper.trenchTable(result.trenches, true);
  }
});

document.getElementById('commit').addEventListener('click', async () => {
  const note = document.getElementById('note').value;
  const result = await api('/api/commit', { poles: stepper.poles(), note, events: stepper.declared }, 'Stepping 10 Myr…');
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
