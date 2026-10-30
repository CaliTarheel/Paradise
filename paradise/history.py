"""A world's history as a log of decisions, stepped a fixed interval at a time.

A project is the painted start and a list of steps. Each step is the stage pole of every plate for
the next `step_myr` (and a note saying why), which is all a person decides in this phase. Replaying
the log from the start rebuilds the world exactly, so the log is the history and the world is
derived from it. A step is never rewritten: going back means truncating the log to an earlier step
and continuing from there.

Poles follow the GPlates convention: latitude and longitude of the pole, and a rate in degrees per
Myr, positive anticlockwise seen from above the pole. A plate without a pole for a step stands still.
"""
from __future__ import annotations

import base64
import copy
import hashlib
import json
import math
import os
import pickle
import zlib
from dataclasses import dataclass, field

import numpy as np
import spherely

from paradise import events as declared
from paradise import collide, orogeny, plumes, rifting, terranes
from paradise.core import sphere
from paradise.core.world import INITIATION_KM
from paradise.paint import world_from_grid

STEP_MYR = 10.
# A step is advanced in fixed sub-steps of this length. The core resolves each sub-step's overlaps
# and gaps along stage lines; shorter sub-steps only move slivers between neighbours, at several
# times the cost, and are not worth it for a worldbuilding history. A fast plate therefore moves a
# few hundred km in one sub-step, and accepts the coarser splitting that comes with it.
SUBSTEP_MYR = 2.5
# With rigid continents, colliding continents weld once they have closed this far: about one
# sub-step's contact, so they meet and stop instead of pushing a hundred km into each other.
RIGID_SUTURE_KM = 5.
SPEED_LIMIT_CM_YR = 20.


def _atomic_geographies():
    """Let deepcopy share spherely geographies: they are never changed in place."""
    for name in dir(spherely):
        kind = getattr(spherely, name)
        if isinstance(kind, type) and issubclass(kind, getattr(spherely, 'Geography', ())):
            copy._deepcopy_dispatch[kind] = lambda value, memo: value


_atomic_geographies()


def pole_vector(lat: float, lon: float, rate_deg_myr: float) -> np.ndarray:
    """Rotation vector (radians per Myr) of a stage pole."""
    return sphere.from_lonlat(lon, lat) * math.radians(rate_deg_myr)


def vector_pole(vector) -> tuple:
    """(lat, lon, rate in degrees per Myr) of a rotation vector, with a positive rate."""
    vector = np.asarray(vector, float)
    rate = float(np.linalg.norm(vector))
    if rate == 0.:
        return 90., 0., 0.
    lon, lat = sphere.to_lonlat(vector / rate)
    return float(lat), float(lon), math.degrees(rate)


def speed_cm_yr(vector, point) -> float:
    """Speed at `point` (unit vector) of a plate turning by `vector` (radians per Myr), cm/yr."""
    km_per_myr = float(np.linalg.norm(np.cross(vector, point))) * sphere.RADIUS_KM
    return km_per_myr / 10.   # 1 km/Myr = 0.1 cm/yr


def encode_grid(grid: np.ndarray) -> dict:
    grid = np.ascontiguousarray(grid, np.uint16)
    return dict(shape=list(grid.shape), dtype='uint16',
                data=base64.b64encode(zlib.compress(grid.tobytes(), 9)).decode('ascii'))


def decode_grid(record: dict) -> np.ndarray:
    raw = zlib.decompress(base64.b64decode(record['data']))
    return np.frombuffer(raw, np.dtype(record.get('dtype', 'uint16'))).reshape(record['shape']).astype(np.int64)


@dataclass
class Snapshot:
    """The world at the end of a step, and what happened during it."""
    time_myr: float
    world: object
    cratons: list
    frames: list          # [(time_myr, {plate: 3x3 frame})], one per sub-step, for the rotation file
    entries: list         # the core's ledger entry of each sub-step
    warnings: list = field(default_factory=list)
    trenches: list = field(default_factory=list)   # paradise.events.trenches at the end of the step
    features: list = field(default_factory=list)   # paradise.terranes.Feature: arcs and terranes
    rifts: list = field(default_factory=list)      # paradise.rifting: lines plates were split along
    new_poles: dict = field(default_factory=dict)  # plate made this step -> the pole it was given
    collisions: list = field(default_factory=list)  # paradise.events.collisions at the end of the step
    orogens: list = field(default_factory=list)     # paradise.orogeny.belts at the end of the step
    hotspots: list = field(default_factory=list)    # paradise.plumes: every hotspot so far, with its trail


@dataclass
class Project:
    name: str
    plates: np.ndarray
    crust: np.ndarray
    start_age_ma: float = 2000.
    step_myr: float = STEP_MYR
    ocean_age_myr: float = 60.
    plate_names: dict = field(default_factory=dict)
    # [{'poles': {plate: [lat, lon, rate]}, 'events': [declaration], 'note': str}]
    steps: list = field(default_factory=list)
    # Whether a converging boundary becomes a trench by itself once it has closed INITIATION_KM.
    # Without it, every trench is declared (paradise.events).
    auto_trenches: bool = True
    # Whether a trench starting between a continental and an oceanic plate always has the continental
    # one overriding (Worldbuilding Pasta's rule for its starting plates).
    continents_override: bool = False
    # Whether only two continental plates weld when they collide (Scotese's rule IV).
    weld_continents_only: bool = False
    # Whether continents are rigid, as in a hand-built reconstruction: they weld where they touch
    # instead of shortening into each other, and nothing escapes sideways along their coasts.
    rigid_continents: bool = False
    # How many Myr of ocean floor each ocean crust piece spans (its isochrons are this far apart).
    ocean_band_myr: float = 10.
    # The largest continental piece that accretes as a terrane instead of colliding (0: none do).
    terrane_max_km2: float = terranes.TERRANE_MAX_KM2

    # ------------------------------------------------------------------ files
    def to_json(self) -> dict:
        return dict(format='paradise-project', version=1, name=self.name, start_age_ma=self.start_age_ma,
                    step_myr=self.step_myr, ocean_age_myr=self.ocean_age_myr, auto_trenches=self.auto_trenches,
                    continents_override=self.continents_override, weld_continents_only=self.weld_continents_only,
                    rigid_continents=self.rigid_continents, ocean_band_myr=self.ocean_band_myr,
                    terrane_max_km2=self.terrane_max_km2,
                    plate_names={str(k): v for k, v in self.plate_names.items()},
                    plates=encode_grid(self.plates), crust=encode_grid(self.crust), steps=self.steps)

    @classmethod
    def from_json(cls, data: dict) -> 'Project':
        if data.get('format') != 'paradise-project':
            raise ValueError('not a Paradise project')
        return cls(data['name'], decode_grid(data['plates']), decode_grid(data['crust']),
                   float(data.get('start_age_ma', 2000.)), float(data.get('step_myr', STEP_MYR)),
                   float(data.get('ocean_age_myr', 60.)),
                   {int(k): v for k, v in data.get('plate_names', {}).items()}, list(data.get('steps', [])),
                   bool(data.get('auto_trenches', True)),
                   bool(data.get('continents_override', False)),
                   bool(data.get('weld_continents_only', False)),
                   bool(data.get('rigid_continents', False)),
                   float(data.get('ocean_band_myr', 10.)),
                   float(data.get('terrane_max_km2', terranes.TERRANE_MAX_KM2)))

    def save(self, path):
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_json(), f)

    @classmethod
    def load(cls, path) -> 'Project':
        with open(path, encoding='utf-8') as f:
            return cls.from_json(json.load(f))

    def age_ma(self, time_myr: float) -> float:
        """Geological age of model time `time_myr` (0 at the start)."""
        return self.start_age_ma - time_myr


CACHE_VERSION = 9


class History:
    """A project and the world after each of its steps, kept so any step can be shown or gone back to.

    Replaying a long history takes a while, so the snapshots can be kept in a cache file beside the
    project (`cache_path`). The cache is only a shortcut: each cached snapshot is used only if it was
    made from the same start and the same steps, and anything else is replayed.
    """

    def __init__(self, project: Project, cache_path=None):
        self.project = project
        self.cache_path = cache_path
        cached = self._read_cache()
        if cached and cached[0][0] == self._key(0):
            self.snapshots = [cached[0][1]]
        else:
            world, cratons = world_from_grid(project.plates, project.crust, project.ocean_age_myr,
                                             project.plate_names)
            self.snapshots = [Snapshot(0., world, cratons, [(0., _frames(world))], [])]
        for k, step in enumerate(project.steps, start=1):
            if k < len(cached) and cached[k][0] == self._key(k) and len(self.snapshots) == k:
                self.snapshots.append(cached[k][1])
            else:
                self.snapshots.append(self._advance(self.snapshots[-1], step['poles'], step.get('events', [])))

    def _key(self, k: int) -> str:
        """What snapshot `k` depends on: the painting, the settings and the first `k` steps."""
        p = self.project
        digest = hashlib.sha256()
        digest.update(json.dumps([CACHE_VERSION, SUBSTEP_MYR, p.step_myr, p.ocean_age_myr, p.auto_trenches, p.continents_override, p.weld_continents_only, p.rigid_continents, p.ocean_band_myr,
                                  p.terrane_max_km2,
                                  list(p.plates.shape)]).encode())
        digest.update(np.ascontiguousarray(p.plates, np.int64).tobytes())
        digest.update(np.ascontiguousarray(p.crust, np.int64).tobytes())
        digest.update(json.dumps([[s['poles'], s.get('events', [])] for s in p.steps[:k]], sort_keys=True).encode())
        return digest.hexdigest()

    def _read_cache(self) -> list:
        if not self.cache_path or not os.path.exists(self.cache_path):
            return []
        try:
            with open(self.cache_path, 'rb') as f:
                return pickle.load(f)
        except Exception:  # noqa: BLE001 - a cache that cannot be read is replayed instead
            return []

    def save_cache(self):
        if not self.cache_path:
            return
        entries = [(self._key(k), snapshot) for k, snapshot in enumerate(self.snapshots)]
        temporary = self.cache_path + '.part'
        with open(temporary, 'wb') as f:
            pickle.dump(entries, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, self.cache_path)

    @property
    def current(self) -> Snapshot:
        return self.snapshots[-1]

    def preview(self, poles: dict, events: list = ()) -> Snapshot:
        """The next step with these poles and declarations, without adding it to the history."""
        return self._advance(self.current, poles, [declared.validate(e) for e in events])

    def commit(self, poles: dict, note: str = '', events: list = ()) -> Snapshot:
        poles = {int(p): [float(v) for v in pole] for p, pole in poles.items()}
        events = [declared.validate(e) for e in events]
        snapshot = self._advance(self.current, poles, events)
        step = dict(poles={str(p): v for p, v in poles.items()}, note=note)
        if events:
            step['events'] = events
        self.project.steps.append(step)
        self.snapshots.append(snapshot)
        return snapshot

    def back_to(self, step: int):
        """Drop every step after `step` (0 is the painted start)."""
        if not 0 <= step < len(self.snapshots):
            raise ValueError('no step %d' % step)
        del self.snapshots[step + 1:]
        del self.project.steps[step:]

    # ------------------------------------------------------------------ stepping
    def _advance(self, before: Snapshot, poles: dict, events: list = ()) -> Snapshot:
        world = copy.deepcopy(before.world)
        world.initiation_km = INITIATION_KM if self.project.auto_trenches else math.inf
        world.terrane_max_km2 = self.project.terrane_max_km2
        world.continents_override = self.project.continents_override
        world.weld_continents_only = self.project.weld_continents_only
        world.age_band_myr = self.project.ocean_band_myr
        if self.project.rigid_continents:
            world.suture_km, world.escape_fraction = RIGID_SUTURE_KM, 0.
        cratons = list(before.cratons)
        features = copy.deepcopy(before.features)
        vectors = {int(p): pole_vector(*pole) for p, pole in poles.items()}
        frames, entries = [], []
        # Declarations apply to the world as it was shown when they were made, before anything moves.
        notes = declared.apply(world, events)
        rifts = copy.deepcopy(before.rifts)
        split, new_poles = [], {}
        for event in events:
            if event['kind'] not in declared.SPLITS:
                continue
            try:
                new, cratons, record = _split(world, event, vectors, cratons, features, rifts)
            except ValueError as e:
                notes.append(dict(kind='split_refused', text='The declared %s was not made: %s.'
                                  % (event['kind'].replace('backarc', 'back-arc'), e)))
                continue
            if 'pole' in event:
                vectors[new] = pole_vector(*event['pole'])
                new_poles[new] = list(event['pole'])
            split.append(record)
        hotspots = copy.deepcopy(before.hotspots)
        for event in events:
            if event['kind'] == 'hotspot':
                split.append(plumes.place_hotspot(hotspots, event, world.time_myr))
            elif event['kind'] == 'end_hotspot':
                found = plumes.end_hotspot(hotspots, event['hotspot'], world.time_myr)
                notes += [f for f in found if 'text' in f]
                split += [f for f in found if 'text' not in f]
            elif event['kind'] == 'lip':
                split.append(plumes.emplace_lip(world, features, event, world.time_myr))
                if event.get('with_hotspot', True):
                    split.append(plumes.place_hotspot(hotspots, dict(point=event['point'],
                                                                     name=event.get('name')), world.time_myr))
        if any(e['kind'] == 'lip' for e in events):
            world.regions = world._from_crust()
        # The new poles take effect now: type every boundary from them before the first sub-step, so
        # its opening is booked as ridge (not as a sliver of a boundary typed by the old motion).
        # Its events (a declared weld happens here) are the step's first.
        split += world._evolve(_rotations(world, vectors), 0., dict(continent_arrived=[]))
        remaining = self.project.step_myr
        while remaining > 1e-9:
            rotations = _rotations(world, vectors)
            sub = min(remaining, SUBSTEP_MYR)
            # Rigid continents about to touch weld where they are, before the core can overlap them.
            welded = collide.weld_on_approach(world, rotations, sub) if self.project.rigid_continents else []
            if welded:
                rotations = _rotations(world, vectors)
            # Terranes about to reach a trench join the overriding plate before the core can collide them.
            accreted = terranes.accrete(world, features, rotations, sub)
            entry = world.step(rotations, sub)
            entry['events'] = welded + accreted + entry['events']
            remaining -= sub
            features = terranes.tidy(world, features)
            plumes.rehome_marks(world, hotspots)
            plumes.mark_trails(world, hotspots, world.time_myr)
            entries.append(_compact(entry))
            frames.append((world.time_myr, _frames(world)))
            cratons = _rehome(cratons, world)
        moving = _rotations(world, vectors)
        trenches = declared.trenches(world, moving, self.project.step_myr)
        born = terranes.grow_arcs(world, trenches, before.trenches, features, world.time_myr,
                                 self.project.step_myr)
        if born:
            entries[-1]['events'] = entries[-1].get('events', []) + born
        snapshot = Snapshot(world.time_myr, world, cratons, frames, entries)
        snapshot.features = features
        snapshot.rifts = rifts
        snapshot.new_poles = new_poles
        if split:
            entries[0]['events'] = _compact(dict(events=split))['events'] + entries[0].get('events', [])
        snapshot.trenches = trenches
        snapshot.collisions = declared.collisions(world)
        snapshot.orogens = orogeny.belts(world)
        snapshot.hotspots = hotspots
        kept = [dict(kind='kept_apart', plates=c['plates'], text=(
            'Plates %d and %d collide along %d km but are kept apart; the continent between them goes '
            'on shortening and thickening.' % (c['plates'][0], c['plates'][1], c['length_km'])))
            for c in snapshot.collisions if c['kept_apart']]
        snapshot.warnings = (notes + review(world, moving) + declared.trench_warnings(snapshot.trenches, moving)
                             + kept)
        return snapshot


def _rotations(world, vectors: dict) -> dict:
    """Every plate's rotation: its own pole, or the pole of the plate it came from (a rift without a
    pole of its own, a piece coming apart during the step), or standing still."""
    rotations = {}
    for p in world.plates():
        ancestor = p
        while ancestor not in vectors and ancestor in world.parents:
            ancestor = world.parents[ancestor]
        rotations[p] = vectors.get(ancestor, np.zeros(3))
    return rotations


def _split(world, event, vectors, cratons, features, rifts):
    """Apply a rift or back-arc declaration; returns (new plate, cratons, event record)."""
    full = {p: vectors.get(p, np.zeros(3)) for p in world.plates()}
    if event['kind'] == 'rift':
        plate = event['plate']
        piece = rifting.rift_piece(world, plate, event['points'], event['new_side'])
        line = sphere.from_lonlat(*np.asarray(event['points'], float).T)
        return rifting.split_off(world, plate, piece, full, 'Plate %d rift' % plate, cratons, features, rifts, line)
    if event['kind'] == 'detach':
        plate = event['plate']
        piece = rifting.old_ocean_piece(world, plate, event['min_age_myr'])
        return rifting.split_off(world, plate, piece, full, 'Plate %d old ocean' % plate, cratons, features, rifts)
    over = event['overriding']
    piece = rifting.backarc_piece(world, event['points'], event['plates'], over)
    return rifting.split_off(world, over, piece, full, 'Plate %d forearc' % over, cratons, features, rifts)


def _frames(world) -> dict:
    return {p: world.frame(p).copy() for p in world.plates()}


def _compact(entry: dict) -> dict:
    """The parts of a core ledger entry worth keeping in the history (no geometry)."""
    keep = {}
    for key, value in entry.items():
        if isinstance(value, (int, float, str)) or key in ('gained', 'lost'):
            keep[key] = value
        elif key == 'events':
            keep[key] = [{k: (v if isinstance(v, (int, float, str, type(None), list)) else
                              list(v) if isinstance(v, tuple) else str(v)) for k, v in e.items()} for e in value]
    return keep


def _rehome(cratons: list, world) -> list:
    """Cratons whose plate has gone move to the plate now under them, in that plate's frame."""
    plates = set(world.plates())
    result = []
    for plate, shape in cratons:
        if plate in plates:
            result.append((plate, shape))
            continue
        # The craton's last frame is unknown once its plate has gone, so its plate's last frame was
        # taken over by whatever absorbed it: look for the plate whose region holds the craton.
        frame = world.frames.get(plate, np.eye(3))
        centre = sphere.rotate(_centre(shape), frame)
        lon, lat = sphere.to_lonlat(centre)
        point = spherely.create_point(float(lon), float(lat))
        new = next((p for p in plates if spherely.contains(world.regions[p], point)), None)
        if new is None:
            continue
        moved = sphere.rotate_geography(shape, world.frame(new).T @ frame)
        result.append((new, moved))
    return result


def _centre(shape) -> np.ndarray:
    c = spherely.centroid(shape)
    return sphere.from_lonlat(spherely.get_x(c), spherely.get_y(c))


def review(world, vectors: dict) -> list:
    """Warnings about a step: plates faster than Earth's fastest, and the core's own checks."""
    warnings = []
    for p in world.plates():
        v = vectors.get(p)
        if v is None or not np.any(v):
            continue
        rings = sphere.polygons(world.regions[p])
        if not rings:
            continue
        points = np.vstack([r for poly in rings for r in poly])
        fastest = max(speed_cm_yr(v, q) for q in points)
        if fastest > SPEED_LIMIT_CM_YR:
            warnings.append(dict(kind='speed', plate=p,
                                 text='Plate %d moves up to %.1f cm/yr; Earth\'s fastest plates move about %d.'
                                      % (p, fastest, SPEED_LIMIT_CM_YR)))
    return warnings
