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

Each step reports its trenches: length, closing rate, the age of the slab, and the area consumed
per step. The warnings come from Scotese's rules of thumb:
- a plate pulling away from its trench;
- continent arriving at a trench;
- a long trench whose downgoing plate barely moves, since plates carried by a slab move fast.

## The painted start

The editor paints two equirectangular grids at 0.5 degrees: plate numbers, and crust (ocean,
continent, craton). `paradise/paint.py` traces each plate-and-crust piece along cell edges into rings
with the piece on their left, so a plate larger than a hemisphere is still the region it encloses,
then simplifies all the rings together so that neighbours still share edges and the plates still
tile the sphere. Cratons are continent in the core and are kept as separate outlines in their
plate's frame.

## GPlates export

Age is the project's start age minus model time. A plate's total rotation at an age is its frame
then, relative to the anchor plate 000; the rotation file has a row for every plate at every
sub-step. Every feature is a snapshot of one step, valid until the next step, with its geometry
carried back into its plate's frame, so GPlates moves it smoothly between steps. The layers follow
Artifexia's: Plates, Continents, Cratons, Ocean crust (named by age), Divergent, Convergent and
Transform.

GPlates is the viewer and the place for final polish. Nothing edited there comes back into Paradise.

## Plan

1. **Done:** paint a start, step it with poles, go back, export to GPlates.
2. Isochrons as lines, the plate hierarchy in the rotation file, and branching to a new project.
3. Events a person declares.
   - **Done:** start subduction with a chosen polarity, flip it, end it.
   - **To do:** rift along a drawn line, trench jump, weld or keep apart, plume or LIP, calve a
     microcontinent, back-arc spreading.
4. An advisor that proposes poles and events from the Scotese rules of thumb, and an autopilot that
   runs several steps and stops at any collision, rift, plate loss or warning.
5. Secondary features: hotspot trails, LIPs, island arcs, back-arcs, accreted terranes, and
   orogenies that age from active to former to old.
6. Elevation and maps, an optional goSPL erosion pass, and Aerie's 72 by 46 ROCKE-3D input.
