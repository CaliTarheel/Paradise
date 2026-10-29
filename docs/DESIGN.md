# Paradise design

Paradise writes a plausible plate-tectonic history for a fictional world. It is not a physics
solver. A person chooses how the plates move, a fixed step at a time; the engine does the
bookkeeping that hand-built GPlates reconstructions (such as the Artifexia project) do by hand, and
exports the result to GPlates.

## Why not a solver

The earlier engines (Deep-Time, Lite, Super-Light-Time, Poly-Time's dynamics, World Crucible) let
forces choose the motions and events. Most of their effort went into collision and continental
entry, slab lifecycles and geometry repair, and the histories still came out unlike Earth's (for
example, mean plate speeds near 1 cm/yr and mean ocean ages near 150 Myr). Choosing the motions
directly removes the forces, and with them most of that work.

## History is a log of decisions

A project is the painted start plus a list of steps. Each step holds one stage pole per plate for
the next 10 Myr and a note saying why. Replaying the log rebuilds the world exactly, so the log is
the history; the world is derived from it. A step is never edited: going back truncates the log, and
the history continues from there.

Within a step every plate turns about one pole, so a step maps directly onto rotation-file rows.

## Steps and sub-steps

The decision step is 10 Myr. The engine advances it in fixed sub-steps of 2.5 Myr. Shorter
sub-steps only move slivers between neighbouring plates at several times the cost, and are not
worth it for a worldbuilding history; a fast plate therefore moves a few hundred kilometres in one
sub-step, with the coarser overlap and gap splitting that implies.

## The core

`paradise/core` is Poly-Time's regions-first plate model (see its `PROVENANCE.md`). Plates are
their crust, stored in each plate's own frame as GPlates stores features. Each sub-step carries every
plate by its rotation, resolves overlaps and gaps by boundary type (trenches consume, ridges create
ocean, everything else splits along the half-stage line), and derives the boundaries again, keeping
each segment's state by position. The crust ledger is checked every sub-step: area enters only as
new ocean and leaves only at boundaries.

Boundary types follow the motion. A converging boundary becomes a trench once it has closed 100 km,
with ocean going under continent and older ocean under younger; a trench a continent reaches becomes
a collision, and a collision that has closed far enough welds the two plates.

## Declared events

A step can carry declarations as well as poles. Each declaration names two plates and a stretch of
the boundary between them, as points picked on the map at the start of the step. It applies to
every segment of that boundary within 60 km of the stretch. Replaying the log reaches the same world
at the start of the step, so the same segments are found again.

Declarations apply before anything moves, and before the boundaries are re-typed from the step's
poles. The core never re-types a trench from the motion, and it carries each segment's type and
overriding plate as the boundary moves. A declared trench therefore lasts until it is ended, or
until a continent reaches it and it becomes a collision.

- **Subduction** makes the segments a trench with the named plate overriding. Declaring it again
  with the other plate over it flips the polarity.
- **Ending subduction** makes the segments converging again, with their closing so far forgotten.
  They are then typed from the motion.

A project can turn off the core's own trench initiation, which otherwise makes a converging
boundary a trench once it has closed 100 km. With it off, every trench is declared. A converging
boundary with no trench is then resolved along its half-stage line, and the ocean lost there is
booked as compression.

Three more declarations decide whether colliding plates become one. Left to the core, a collision
that has closed 100 km welds its two plates.
- **Weld** makes every segment between the two plates a collision closed that far, so the core
  welds them in the re-typing before anything moves.
- **Keep apart** adds the pair to the world's `kept_apart` set, which the core's welding skips. It
  lasts until **allow weld** removes it.

Each step reports its trenches: length, closing rate, the age of the slab, and the area consumed
per step. The warnings come from Scotese's rules of thumb:
- a plate pulling away from its trench;
- continent arriving at a trench;
- a long trench whose downgoing plate barely moves, since plates carried by a slab move fast.

## Rifts and back-arcs

`paradise/rifting.py` splits a plate in two at the start of a step. The new plate starts in its
parent's frame, so its crust keeps its coordinates and moves with the parent until it is given a
pole. Its number is one more than any plate the world has had, which replays the same every time.

- **A rift** is a line drawn across the plate. The plate less a 1 km cut along the line falls into
  parts, and the parts on the named side become the new plate. A line that does not cross the plate
  from edge to edge leaves only one side, and is refused with a note.
- **A back-arc** is the band of the overriding plate from just outboard of a trench to 250 km
  inboard. That is past the arc, at 150 km, with room for the basin to open.

The core carries a boundary's state only between the same two plates. So, before the boundaries are
derived again, the parent's boundary segments along the new plate's edge are relabelled to the new
plate. A trench a forearc takes over stays a trench, with the forearc overriding it. Crust,
cratons, arcs and terranes, sutures and earlier rift lines on the new side go with the new plate.

## Island arcs and terranes

`paradise/terranes.py` adds two rules on top of the trenches. Both are bookkeeping; neither is
mechanics.

**Arcs.** At the end of a step, a trench of the same plates and polarity as in the previous step,
closing at least 1 cm/yr, raises a band 30 km either side of a line 150 km inboard on the
overriding plate.
- The band's ocean becomes continental crust 25 km thick: an island arc.
- Its continent is recorded as a continental arc, adding no crust.

Arc crust that the band still covers thickens 0.2 km per Myr, up to 35 km, as an arc matures.

Requiring the trench to have lasted a step stands in for the few million years a slab takes to
reach the depth where arc magmatism starts. The arc is on the overriding plate, so it stays behind
when its trench moves away (Scotese: arcs do not ride their trenches across oceans).

**Terranes.** Before each sub-step, any continental piece of at most the project's terrane size
(1 million km² unless changed) that will reach,
within that sub-step's motion plus 50 km, a trench its plate goes down is given to the overriding
plate. The ocean around it, out to the trench, goes with it.
- The handover happens before the core moves anything. Otherwise the core would turn the trench
  into a collision as soon as 500 km² of continent arrived, and later weld the two plates.
- The boundaries are derived again, keeping their state by position. The boundary round the
  terrane can bulge further outboard than the core carries a segment's state (300 km), so every
  segment between the two plates along the accreted zone is made trench again, with the same
  overriding plate. The trench thus steps outboard of the terrane.

Arcs and terranes are features kept beside the world, each in its plate's frame, as cratons are.
They keep their identity and story while the core's crust pieces merge and split. After each
sub-step a feature is clipped to its plate's continental crust, and moved to the plate under it if
its own plate has gone.

## Orogenies, hotspots and LIPs

**Orogenies** (`paradise/orogeny.py`) are derived from the world at the end of each step, not
carried. Each belt is the continental crust within 250 km of a line:
- **active**: along current collisions;
- **former**: along sutures made less than 100 Myr ago;
- **old**: along older sutures.

The suture lines the core keeps when it welds two plates are what remembers where an orogeny was.

**Plumes** (`paradise/plumes.py`) are declared at a point.
- **A hotspot** is fixed in the anchor frame. After every sub-step the plate over it gets a trail
  mark in that plate's frame, so the trail moves with the plate. Marks on a plate that has gone
  move to the plate now under them, as cratons do.
- **A LIP** is a cap emplaced at once.
  - Its ocean becomes 30 km oceanic plateau, continental crust in the core. It is therefore
    buoyant, and accretes or collides under the terrane rules.
  - Its continent is recorded as a flood-basalt province.
  - It leaves a hotspot at its centre unless declared without one.

A pair of plates can be barred from accretion (`no_accretion`, a set on the world that
`terranes.accrete` skips), so that even a small terrane collides there.

## The painted start

The editor paints two equirectangular grids at 0.5 degrees: plate numbers, and crust (ocean,
continent, craton). `paradise/paint.py` traces each plate-and-crust piece along cell edges into rings
with the piece on their left, so a plate larger than a hemisphere is still the region it encloses,
then simplifies all the rings together so that neighbours still share edges and the plates still
tile the sphere. No edge is left longer than two cells (about 110 km); the core's own simplifying of moved
outlines allows 150 km. A declaration applies to whole segments, so that is how closely one can be
placed. The finer edges cost nothing measurable
per step. Cratons are continent in the core and are kept as separate outlines in their
plate's frame.

## GPlates export

Age is the project's start age minus model time. A plate's total rotation at an age is its frame
then, and the rotation file has a row for every plate at every sub-step.
- A plate made by a rift or a back-arc rotates relative to its parent, as in a hand-built tree,
  when the parent exists at every one of its rows.
- Every other plate rotates relative to the anchor plate 000. Every feature is a snapshot of one step, valid until the next step, with its geometry
carried back into its plate's frame, so GPlates moves it smoothly between steps. The layers follow
Artifexia's: Plates, Continents, Cratons, Ocean crust (named by age), Divergent, Convergent and
Transform.

GPlates is the viewer and the place for final polish. Nothing edited there comes back into Paradise.

## Plan

1. **Done:** paint a start, step it with poles, go back, export to GPlates.
2. Isochrons as lines, the plate hierarchy in the rotation file, and branching to a new project.
3. Events a person declares.
   - **Done:** start subduction with a chosen polarity, flip it, end it.
   - **Done, as rules rather than declarations:** island and continental arcs behind lasting
     trenches; terranes that accrete instead of colliding.
   - **Done:** rift along a drawn line; open a back-arc behind a trench.
   - **Done:** weld or keep apart; hotspots and LIPs; accretion barred per pair.
   - **To do:** trench jump.
4. An advisor that proposes poles and events from the Scotese rules of thumb, and an autopilot that
   runs several steps and stops at any collision, rift, plate loss or warning.
5. **Done:** secondary features: hotspot trails, LIPs, island arcs, back-arcs, accreted terranes,
   and orogenies that age from active to former to old.
6. Elevation and maps, an optional goSPL erosion pass, and Aerie's 72 by 46 ROCKE-3D input.
