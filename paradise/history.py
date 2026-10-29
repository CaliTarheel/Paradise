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
from paradise.core import sphere
from paradise.core.world import INITIATION_KM
from paradise.paint import world_from_grid

STEP_MYR = 10.
# A step is advanced in fixed sub-steps of this length. The core resolves each sub-step's overlaps
# and gaps along stage lines; shorter sub-steps only move slivers between neighbours, at several
# times the cost, and are not worth it for a worldbuilding history. A fast plate therefore moves a
# few hundred km in one sub-step, and accepts the coarser splitting that comes with it.
SUBSTEP_MYR = 2.5
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

    # ------------------------------------------------------------------ files
    def to_json(self) -> dict:
        return dict(format='paradise-project', version=1, name=self.name, start_age_ma=self.start_age_ma,
                    step_myr=self.step_myr, ocean_age_myr=self.ocean_age_myr, auto_trenches=self.auto_trenches,
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
                   bool(data.get('auto_trenches', True)))

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


CACHE_VERSION = 2


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
        digest.update(json.dumps([CACHE_VERSION, SUBSTEP_MYR, p.step_myr, p.ocean_age_myr, p.auto_trenches,
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
        cratons = list(before.cratons)
        vectors = {int(p): pole_vector(*pole) for p, pole in poles.items()}
        frames, entries = [], []
        # Declarations apply to the world as it was shown when they were made, before anything moves.
        notes = declared.apply(world, events)
        # The new poles take effect now: type every boundary from them before the first sub-step, so
        # its opening is booked as ridge (not as a sliver of a boundary typed by the old motion).
        world._evolve({p: vectors.get(p, np.zeros(3)) for p in world.plates()}, 0., dict(continent_arrived=[]))
        remaining = self.project.step_myr
        while remaining > 1e-9:
            rotations = {p: vectors.get(p, np.zeros(3)) for p in world.plates()}
            # A plate made during the step (a piece coming apart) moves with the plate it came from.
            for p in world.plates():
                ancestor = p
                while ancestor not in vectors and ancestor in world.parents:
                    ancestor = world.parents[ancestor]
                rotations[p] = vectors.get(ancestor, np.zeros(3))
            sub = min(remaining, SUBSTEP_MYR)
            entry = world.step(rotations, sub)
            remaining -= sub
            entries.append(_compact(entry))
            frames.append((world.time_myr, _frames(world)))
            cratons = _rehome(cratons, world)
        snapshot = Snapshot(world.time_myr, world, cratons, frames, entries)
        snapshot.trenches = declared.trenches(world, vectors, self.project.step_myr)
        snapshot.warnings = notes + review(world, vectors) + declared.trench_warnings(snapshot.trenches, vectors)
        return snapshot


def _frames(world) -> dict:
    return {p: world.frame(p).copy() for p in world.plates()}


def _compact(entry: dict) -> dict:
    """The parts of a core ledger entry worth keeping in the history (no geometry)."""
    keep = {}
    for key, value in entry.items():
        if isinstance(value, (int, float, str)) or key in ('gained', 'lost'):
            keep[key] = value
        elif key == 'events':
            keep[key] = [{k: (v if isinstance(v, (int, float, str, type(None))) else
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
        point = spherely.point(float(lon), float(lat))
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
