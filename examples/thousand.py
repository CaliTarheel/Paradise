"""A thousand-million-year history, by autopilot: the advisor's proposals, committed step by step.

    python examples/thousand.py worlds/thousand.paradise.json --steps 100

It builds a Pangaea start, as the hand-built Artifexia project has (paradise.starts: one
supercontinent on one plate in one ocean plate), with rigid continents, or carries on from the project file
if it exists, and saves the project after every step. A step that fails is tried again without its
declarations, then with every plate at half speed, and the run stops if that fails too, with the
error, so nothing is skipped silently.
"""
import argparse
import os
import resource
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from paradise.advisor import Advisor                      # noqa: E402
from paradise.core import sphere                          # noqa: E402
from paradise.history import History, Project             # noqa: E402
from paradise.starts import pangaea_world                 # noqa: E402


def summary(snapshot):
    world = snapshot.world
    continent = sum(c.area_km2() for c in world.crust if c.kind == 'continent') / sphere.SPHERE_KM2
    ocean = [(c.area_km2(), world.time_myr - c.born_myr) for c in world.crust if c.kind == 'ocean' and c.born_myr is not None]
    age = sum(a * t for a, t in ocean) / max(sum(a for a, _ in ocean), 1.)
    kinds = {}
    for f in snapshot.features:
        kinds[f.kind] = kinds.get(f.kind, 0) + 1
    vertices = sum(len(r) for c in world.crust for poly in sphere.polygons(c.shape) for r in poly)
    return ('plates %d, continent %.0f%%, mean ocean age %.0f Myr, trenches %d, collisions %d, %s, vertices %d'
            % (len(world.plates()), 100 * continent, age, len(snapshot.trenches), len(snapshot.collisions),
               ', '.join('%s %d' % kv for kv in sorted(kinds.items())) or 'no arcs', vertices))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('project')
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()
    if os.path.exists(args.project):
        project = Project.load(args.project)
    else:
        plates, crust = pangaea_world(seed=args.seed)
        project = Project('Thousand', plates, crust, start_age_ma=1000., auto_trenches=True,
                          continents_override=True, weld_continents_only=True, rigid_continents=True,
                          ocean_band_myr=50.)
    t = time.time()
    history = History(project)
    print('ready at %.0f Ma in %.0f s: %s' % (project.age_ma(history.current.time_myr), time.time() - t,
                                             summary(history.current)), flush=True)
    advisor = Advisor(args.seed, rift_at_start=True)
    while len(project.steps) < args.steps:
        poles, events, reasons = advisor.propose(history)
        t = time.time()
        snapshot = None
        for attempt, (p, e) in enumerate(((poles, events), (poles, []),
                                          ({k: [v[0], v[1], v[2] / 2.] for k, v in poles.items()}, []))):
            try:
                snapshot = history.commit(p, 'autopilot: ' + '; '.join(reasons[:3]) if attempt == 0 else
                                          'autopilot, try %d' % (attempt + 1), events=e)
                break
            except Exception:  # noqa: BLE001 - reported, then tried again more gently
                print('step %d attempt %d failed:\n%s' % (len(project.steps) + 1, attempt + 1, traceback.format_exc()),
                      flush=True)
        if snapshot is None:
            print('stopping: three attempts failed', flush=True)
            break
        advisor.committed(snapshot)
        project.save(args.project)
        if len(project.steps) % 10 == 0:
            history.cache_path = args.project + '.cache'
            history.save_cache()   # so maps and the export need not replay the whole run
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
        print('%4.0f Ma  step %3d  %5.1f s  %s%s  (%.0f MB)' % (
            project.age_ma(snapshot.time_myr), len(project.steps), time.time() - t, summary(snapshot),
            '  [%s]' % events[0]['kind'] if events else '', rss), flush=True)


if __name__ == '__main__':
    main()
