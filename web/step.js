// Stepping the history: a stage pole per plate, preview, commit, go back, export.
'use strict';

const KIND_COLOURS = {
  ridge: '#d0362b', trench: '#2a5fae', transform: '#2f8a4f', converging: '#c78a12', collision: '#7a4a1c',
};

const OROGEN_COLOURS = { active: '#c0392b', former: '#8e5a2a', old: '#9c9c8a' };
const FEATURE_COLOURS = { island_arc: '#e2561f', continental_arc: '#b02a8a', terrane: '#6b3fa0', lip: '#7d1d6b' };

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
  // A rift line being drawn: lon/lat points.
  rift: [],
  // A plume being placed: [lon, lat] or null.
  plume: null,

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
    this.rift = [];
    this.riftPanel();
    this.plume = null;
    this.plumePanel();
    this.hotspotList(state.lines.hotspots);
    this.table();
    this.history();
    this.notes(state.warnings, state.events);
    this.trenchTable(state.trenches);
    this.featureList(state.features);
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

  hotspotList(hotspots) {
    const list = document.getElementById('hotspots');
    list.innerHTML = '';
    for (const h of hotspots || []) {
      const li = document.createElement('li');
      li.textContent = `${h.name} (${h.point[0].toFixed(1)}, ${h.point[1].toFixed(1)})${h.active ? '' : ', ended'}`;
      if (h.active) {
        const b = document.createElement('button');
        b.textContent = 'End';
        b.addEventListener('click', () => {
          this.declared.push({ kind: 'end_hotspot', hotspot: h.ident, name: h.name });
          this.preview = null; this.declaredList(); this.draw();
        });
        li.appendChild(b);
      }
      list.appendChild(li);
    }
    if (!list.children.length) list.innerHTML = '<li class="hint">None. Place a plume on the map to add one.</li>';
  },

  plumePanel() {
    const panel = document.getElementById('plume-panel');
    panel.hidden = !this.plume;
    if (!this.plume) return;
    document.getElementById('plume-text').textContent =
      `Plume at ${this.plume[0].toFixed(1)}, ${this.plume[1].toFixed(1)}. A hotspot stays fixed in the mantle and leaves a trail on the plates that pass over it; a LIP erupts at once, making oceanic plateau over ocean and flood basalt over continent.`;
  },

  declarePlume(kind, withHotspot) {
    const event = { kind, point: this.plume.slice() };
    const name = document.getElementById('plume-name').value.trim();
    if (name) event.name = name;
    if (kind === 'lip') { event.radius_km = Number(document.getElementById('plume-radius').value); event.with_hotspot = withHotspot; }
    this.declared.push(event);
    this.plume = null;
    this.preview = null;
    this.plumePanel(); this.declaredList(); this.draw();
  },

  featureList(features) {
    const list = document.getElementById('features');
    list.innerHTML = '';
    for (const f of features || []) {
      const li = document.createElement('li');
      const name = { island_arc: 'Island arc', continental_arc: 'Continental arc', terrane: 'Terrane', lip: 'LIP' }[f.kind] || f.kind;
      li.textContent = `${name} ${f.ident} on plate ${f.plate}, ${Math.round(f.area_km2 / 1000).toLocaleString()} thousand km²: ${f.story}.`;
      list.appendChild(li);
    }
    if (!list.children.length) list.innerHTML = '<li class="hint">None yet. Arcs rise behind trenches that have lasted a step.</li>';
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
    if (!this.pick) { panel.hidden = true; this.newPolePanel(); return; }
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
    const backarc = document.getElementById('open-backarc');
    backarc.hidden = run.kind !== 'trench';
    backarc.onclick = () => this.declare('backarc', run.overriding);
    const kept = (this.state.kept_apart || []).some(([x, y]) => x === a && y === b);
    document.getElementById('keep-apart').hidden = kept;
    document.getElementById('allow-weld').hidden = !kept;
    document.getElementById('weld').onclick = () => this.declarePair('weld');
    document.getElementById('keep-apart').onclick = () => this.declarePair('keep_apart');
    document.getElementById('allow-weld').onclick = () => this.declarePair('allow_weld');
    const barred = (this.state.no_accretion || []).some(([x, y]) => x === a && y === b);
    document.getElementById('no-accretion').hidden = barred || run.kind !== 'trench';
    document.getElementById('allow-accretion').hidden = !barred;
    document.getElementById('no-accretion').onclick = () => this.declarePair('no_accretion');
    document.getElementById('allow-accretion').onclick = () => this.declarePair('allow_accretion');
    this.newPolePanel();
    document.getElementById('pick-clear').onclick = () => { this.pick = null; this.pickPanel(); this.draw(); };
  },

  newPole() {
    const v = ['lat', 'lon', 'rate'].map((k) => Number(document.getElementById(`new-pole-${k}`).value));
    return v[2] ? v : null;
  },

  newPolePanel() {
    const trench = this.pick && this.pick.run.kind === 'trench';
    document.getElementById('new-pole').hidden = !(trench || this.rift.length);
  },

  riftPanel() {
    const panel = document.getElementById('rift-panel');
    panel.hidden = !this.rift.length;
    this.newPolePanel();
    if (!this.rift.length) return;
    const select = document.getElementById('rift-plate');
    const current = Number(select.value);
    select.innerHTML = '';
    for (const p of this.state.plates) {
      const o = document.createElement('option'); o.value = p.id; o.textContent = `Plate ${p.id}`; select.appendChild(o);
    }
    // The plate under the middle of the line, unless one was chosen already.
    const mid = this.rift[Math.floor(this.rift.length / 2)];
    const width = this.state.map_width, height = width / 2;
    const under = this.plateOf[Math.min(height - 1, Math.floor((90 - mid[1]) / 180 * height)) * width
      + Math.min(width - 1, Math.floor((mid[0] + 180) / 360 * width))];
    select.value = this.riftPlateChosen && this.state.plates.some((p) => p.id === current) ? current : under;
    document.getElementById('rift-text').textContent = this.rift.length < 2
      ? 'Click more points. The line must cross the plate from beyond one edge to beyond another.'
      : `Rift line of ${this.rift.length} points, ${Math.round(lengthKm(this.rift)).toLocaleString()} km.`;
    document.getElementById('rift-declare').disabled = this.rift.length < 2;
  },

  declareRift() {
    const side = document.querySelector('input[name=rift-side]:checked').value;
    const event = { kind: 'rift', plate: Number(document.getElementById('rift-plate').value), points: this.rift.slice(), new_side: side };
    const pole = this.newPole();
    if (pole) event.pole = pole;
    this.declared.push(event);
    this.rift = [];
    this.riftPlateChosen = false;
    this.preview = null;
    this.riftPanel(); this.declaredList(); this.trenchTable(this.state.trenches); this.draw();
  },

  declarePair(kind) {
    this.declared.push({ kind, plates: this.pick.run.plates.slice() });
    this.pick = null;
    this.preview = null;
    this.pickPanel(); this.declaredList(); this.trenchTable(this.state.trenches); this.draw();
  },

  declare(kind, overriding) {
    const event = { kind, plates: this.pick.run.plates.slice(), points: this.pickedPoints() };
    if (kind === 'subduction' || kind === 'backarc') event.overriding = overriding;
    if (kind === 'backarc') { const pole = this.newPole(); if (pole) event.pole = pole; }
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
      li.textContent = describeDeclaration(e) + (e.points ? ` (${Math.round(lengthKm(e.points)).toLocaleString()} km)` : '');
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
    for (const o of lines.orogens || []) polyline(ctx, o.points, w, h, OROGEN_COLOURS[o.kind] || '#000', 2, [2, 3]);
    for (const f of lines.features || []) polyline(ctx, f.points, w, h, FEATURE_COLOURS[f.kind] || '#000', 1.5);
    for (const t of lines.trails || []) {
      polyline(ctx, t.points, w, h, 'rgba(0,0,0,0.5)', 4);
      polyline(ctx, t.points, w, h, '#ffd23f', 2);
    }
    for (const s of lines.hotspots || []) {
      const [x, y] = toCanvas(s.point, w, h);
      ctx.beginPath(); ctx.arc(x, y, 6, 0, 2 * Math.PI);
      ctx.fillStyle = s.active ? '#ffd23f' : 'rgba(255,210,63,0.4)'; ctx.fill();
      ctx.strokeStyle = '#000'; ctx.lineWidth = 1.5; ctx.stroke();
    }
    const plumes = this.declared.filter((e) => e.kind === 'hotspot' || e.kind === 'lip').map((e) => e.point);
    if (this.plume) plumes.push(this.plume);
    for (const pt of plumes) {
      const [x, y] = toCanvas(pt, w, h);
      ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.arc(x, y, 9, 0, 2 * Math.PI); ctx.stroke();
    }
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
    // A rift line being drawn, then declared ones.
    const rifts = this.declared.filter((e) => e.kind === 'rift').map((e) => e.points);
    if (this.rift.length) rifts.push(this.rift);
    for (const line of rifts) {
      polyline(ctx, line, w, h, 'rgba(255,255,255,0.95)', 5);
      polyline(ctx, line, w, h, KIND_COLOURS.ridge, 2.5, [6, 4]);
      for (const p of line) { const [x, y] = toCanvas(p, w, h); ctx.fillStyle = '#fff'; ctx.fillRect(x - 3, y - 3, 6, 6); }
    }
    // Declarations waiting for the next step: the stretch, with teeth on the side of the plate that will override.
    for (const e of this.declared.filter((e) => e.kind !== 'rift' && e.points)) {
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
    case 'hotspot_placed': return `Hotspot ${e.name} rose.`;
    case 'hotspot_ended': return `Hotspot ${e.name} died out.`;
    case 'lip_erupted': return `${e.name} erupted: ${Math.round(e.plateau_km2 / 1000)} thousand km² of oceanic plateau and ${Math.round(e.flood_km2 / 1000)} thousand km² of flood basalt.`;
    case 'arc_born': return `An island arc (${Math.round(e.area_km2 / 1000)} thousand km²) rose on plate ${e.plate} behind its trench with plate ${e.plates[0] === e.plate ? e.plates[1] : e.plates[0]}.`;
    case 'accreted': return `A terrane of ${Math.round(e.area_km2 / 1000)} thousand km² from plate ${e.plate} accreted to plate ${e.into}; the trench stepped outboard of it.`;
    case 'rift': return e.name && /forearc/.test(e.name)
      ? `The forearc of plate ${e.plate} split off as plate ${e.new}, to open a back-arc behind it.`
      : `Plate ${e.plate} rifted; plate ${e.new} (${Math.round(e.area_km2 / 1e6 * 10) / 10} million km²) is the new piece.`;
    default: return e.kind ? `${e.kind.replace(/_/g, ' ')}${pair ? ' (' + pair + ')' : ''}` : '';
  }
}

function describeDeclaration(e) {
  const pole = e.pole ? `, new plate pole ${e.pole[0]}, ${e.pole[1]} at ${e.pole[2]}°/Myr` : '';
  if (e.kind === 'hotspot') return `Hotspot${e.name ? ' ' + e.name : ''} at ${e.point[0].toFixed(1)}, ${e.point[1].toFixed(1)}`;
  if (e.kind === 'lip') return `LIP${e.name ? ' ' + e.name : ''} of ${e.radius_km} km at ${e.point[0].toFixed(1)}, ${e.point[1].toFixed(1)}${e.with_hotspot ? ', leaving a hotspot' : ''}`;
  if (e.kind === 'end_hotspot') return `End hotspot ${e.name || e.hotspot}`;
  if (e.kind === 'rift') return `Rift across plate ${e.plate}, new plate on its ${e.new_side}${pole}`;
  const [a, b] = e.plates;
  if (e.kind === 'backarc') return `Back-arc behind the trench between plates ${a} and ${b}, on plate ${e.overriding}${pole}`;
  if (e.kind === 'subduction') return `Trench: plate ${e.overriding === a ? b : a} goes under plate ${e.overriding}`;
  if (e.kind === 'end_subduction') return `End the trench between plates ${a} and ${b}`;
  if (e.kind === 'weld') return `Weld plates ${a} and ${b} into one`;
  if (e.kind === 'keep_apart') return `Keep plates ${a} and ${b} apart when they collide`;
  if (e.kind === 'allow_weld') return `Let plates ${a} and ${b} weld when they collide`;
  if (e.kind === 'no_accretion') return `Continent reaching the trench between plates ${a} and ${b} collides`;
  if (e.kind === 'allow_accretion') return `Let terranes accrete at the trench between plates ${a} and ${b}`;
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
  document.getElementById('map-hint').textContent = {
    boundary: 'Click a boundary to pick it. Click it again to keep only the stretch between your two clicks.',
    rift: 'Click points along the rift, starting and ending beyond the edges of the plate it splits.',
    plume: 'Click where the plume rises.',
    plate: "Click a plate to select its row. Shift-click the map to put the selected plate's pole there.",
  }[mapMode()];
}));

document.getElementById('map').addEventListener('click', (e) => {
  if (!stepper.state) return;
  const canvas = e.target, r = canvas.getBoundingClientRect();
  const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
  if (mapMode() === 'boundary' && !e.shiftKey) {
    stepper.pickAt(fx * canvas.width, fy * canvas.height);
    return;
  }
  if (mapMode() === 'plume' && !e.shiftKey) {
    stepper.plume = [Math.round((fx * 360 - 180) * 100) / 100, Math.round((90 - fy * 180) * 100) / 100];
    stepper.plumePanel(); stepper.draw();
    return;
  }
  if (mapMode() === 'rift' && !e.shiftKey) {
    stepper.rift.push([Math.round((fx * 360 - 180) * 100) / 100, Math.round((90 - fy * 180) * 100) / 100]);
    stepper.riftPanel(); stepper.draw();
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

document.getElementById('rift-declare').addEventListener('click', () => stepper.declareRift());
document.getElementById('rift-undo').addEventListener('click', () => { stepper.rift.pop(); stepper.riftPanel(); stepper.draw(); });
document.getElementById('rift-clear').addEventListener('click', () => {
  stepper.rift = []; stepper.riftPlateChosen = false; stepper.riftPanel(); stepper.draw();
});
document.getElementById('rift-plate').addEventListener('change', () => { stepper.riftPlateChosen = true; });
document.getElementById('plume-hotspot').addEventListener('click', () => stepper.declarePlume('hotspot'));
document.getElementById('plume-lip').addEventListener('click', () => stepper.declarePlume('lip', true));
document.getElementById('plume-lip-only').addEventListener('click', () => stepper.declarePlume('lip', false));
document.getElementById('plume-clear').addEventListener('click', () => { stepper.plume = null; stepper.plumePanel(); stepper.draw(); });
