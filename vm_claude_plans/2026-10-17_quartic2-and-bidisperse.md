**Status:** DRAFT, NOT APPROVED (2026-10-17; revised 2026-10-21 with user decisions, see 'User decisions (2026-10-21)'). User asked to park it at the top of the todo for next session. Nothing executed. Mandatory user code review of every library change before any production run.

# Plan: QUARTIC2.0 validation + production, MTQ3 look, bidisperse (quadratic) tissues

## Context
The experiment needs a regime where the quartic area term dominates. Measured on the new EQ files:
MTQ3 is quartic-dominated (stiffness ratio 6C4(a-1)^2 = 15), QuTs is not (A0 = Abar puts cells at a ~ 1).
In the existing QUARTIC form the quartic coefficient is K*C4, so K = 0 would remove all area elasticity
and divide G's normalization by zero. The user's own QUARTIC2.0 branch,
E_A = 1/2 K A0 [(1-6C4)(a-1)^2 + C4(a-1)^4], has 1/2 f'' = 1 at |a-1| = 1 for every C4, is pure quartic at
C4 = 1/6, and is a double well for C4 > 1/6. The user wants QUARTIC2.0 for three cases, and separately
wants to know whether the model runs with a strongly bidisperse tissue whose EQUILIBRIUM areas are bimodal.

Code survey (Explore agent) found, in `AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp`:
- Forces (AGTissue::ForceFunction -> AGCell::GetPressure) and the Hessian / G code
  (GetPressure, GetNegativeDPressureDArea) already handle QUARTIC2.0 and per-cell A0. No global A0/C4/type there.
- **Bug:** `CalculateFreeEnergy` (:457-500) and `CalculatePressureEnergy` (:898-938) multiply the area sum by
  the LAST cell's `a0` after the loop. Exact for uniform A0, wrong for bidisperse. That energy drives the
  backtracking step test (:644-697) and the equilibrium slope test (:769-825); a mismatch can stall the step
  below 1e-14, which is then written as `_EQ` (:655-660 -> :798-800). The descendant `AGStiffTissue`
  (AGStiffness) does the per-cell sum correctly (`KAM*(A-a0)^2` inside the loop, AGStiffTissue.cpp:48-77).
- **Gap:** per-cell P0 is never read: `(*cptr).P0 = P0;` (:396) uses the header P0. Needed for
  P0_c = p0*sqrt(A0_c). (AGStiffTissue also ignores a P0 column.)
- Apoptosis (MINAREA 0.1, T2 removal) runs when APOPTOSISON = 1. `GetArea` is fabs, so an inverted cell
  reads as positive area.
- Hessian G mode loop returns 0 for floppy modes and silently skips negative eigenvalues.

## User decisions (2026-10-17)
- QUARTIC2.0 cases: A0 = 0.5 Abar & C4 = 1/12; A0 = Abar & C4 = 1/6; A0 = 2 Abar & C4 = 1/2 (double well).
- Double well: **apoptosis OFF**.
- Bidisperse: **quadratic area energy only**; 50/50 mix, mean A0 = Abar; big:small = **1.4:1 and 2:1**;
  **same shape index** P0_c = p0*sqrt(A0_c); the equilibrated areas must come out near-bimodal.

## User decisions (2026-10-21)
- Step 1: mimic `AGStiffTissue::CalculateFreeEnergy` (per-cell a0 inside the loop), unless a correction is warranted.
- QZ sets: K = 1.414214 (the T2ss value) for all three, NOT K proportional to A0.
- Bidisperse source: a quadratic set whose header and geometry match Aq2w. T2ss_OG raw files match it exactly
  (A0 1250, K 1.414214, P0, LAMBDA, ALPHAM; the Vertices section has the same md5 as Aq2w and Qu2w at rng 120), and they are already QUADRATIC.
- Eigenvalue thresholds stay in the library (AGPerimTissue.cpp), never in a Main file. Zero modes are classified with
  an absolute-value bound |lambda| <= threshold; lambda < -threshold is a NEGATIVE mode, reported separately (today it falls
  into the floppy branch).
- Snapshots: borrow the user's MATLAB `MATLAB Code Folder/CreateMovie_HelpperFunctions/PlotSimulationOutput.m`
  (fill() per cell) via a wrapper in the work folder; do not edit repo MATLAB files.

## Step 0 - REMOVED (2026-10-21)
The MTQ3 Hessian and the G vs p0 plot for the A0 = 2 Abar family are already in the C4-redo wrap-up on the todo list. The tissue
snapshots move there too (MATLAB `PlotSimulationOutput.m` wrapper). `plot_tissue.py` is no longer planned.

## Step 1 - Library fix (AGPerimTissue.cpp, one commit)
1. `CalculateFreeEnergy` and `CalculatePressureEnergy`: accumulate `a0 * f(a)` inside the loop and drop the
   trailing `a0` factor (same pattern as AGStiffTissue::CalculateFreeEnergy). Covers every pressure type.
2. Cells reader (:318-405): if the Cells section has a `P0` column, use it per cell; else the header P0 (old
   behavior). Existing `_EQ/_NE` files carry a P0 column equal to the header P0, so restarts are unchanged.
3. Verbose comments on both changes; add a section to `vm_claude_reports/2026-10-15_vertex_model_code_changes_report.md`.
- Regression: on uniform-A0 files (T2ss, QuOw, new QuTs/MTQ3 EQ) the new energy equals `finalEnergy` to ~1e-14
  relative; re-equilibrate 2 QuOw tissues and compare to the 2023 EQ energies (as done for the C4 fix).
- Rebuild the binaries (`agperim_equilibration_claude_chunk`, `agperim_hess_claude`) in a NEW work folder
  `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Quartic2_Bidisperse/` (do not touch the binaries the
  running Quartic_C4fix jobs use).

## Step 2 - QUARTIC2.0 verification (TOP PRIORITY, before any production)
- Extend the independent checks to QUARTIC2.0 and per-cell A0/P0:
  `Quartic_C4fix/analysis/residual_forces/vmcheck.py` (energy, forces) and `Quartic_C4fix/src/RefModel.h` +
  `MainHessFDCheck.cpp` (Hessian, Xi, G vs finite differences).
- Test tissues: MTQ3/Qu2w/Aq2w/Lq2w raw and EQ states converted to QUARTIC2.0 at C4 = 1/12, 1/6, 1/2, plus
  one bidisperse tissue. Pass: library energy = reference; library force = -FD gradient (1e-9 rel);
  Hessian = FD Hessian (~1e-10); G_code = G_FD on non-floppy states.

## Step 3 - QUARTIC2.0 production (3 sets x 112 = 336 tissues)
- Inputs: copy existing raw inputs (same Voronoi geometries, rng 100-160, p0 3.75-3.90) and edit only the
  header C4 / APOPTOSISON and the cell type `QUARTIC` -> `QUARTIC2.0` (Python script, LF endings, exact string):
  | new code | source raw set | A0/Abar | K | C4 | APOPTOSISON |
  |---|---|---|---|---|---|
  | QZ05 | Lq2w | 0.5 | 1.414214 | 1/12 | 1 |
  | QZ10 | Aq2w | 1 | 1.414214 | 1/6 | 1 |
  | QZ20 | Qu2w | 2 | 1.414214 | 1/2 | 0 |
  (2026-10-21: K fixed at the T2ss value for every QZ set, per the user. K/A0 is therefore 2.26e-3 / 1.13e-3 / 5.66e-4,
  not constant across the sets.)
  Editing the raw files is equivalent to re-running MATLAB: every lim-10 set at a given rng has the identical
  Vertices section (agent checked md5), and the generators' K/LAMBDA formulas differ between variants
  (`_SweepC4` ties LAMBDA and ALPHAM to A0), so editing keeps Lambda = 70.71, ALPHAM, box and geometry fixed.
  New codes avoid every name in `PhaseTransitionProject/ExperimentalConditions_CodeNames.csv`; I will add rows there.
  (C4 written with full precision; K edited to 1.414214 in the header.)
- Pilot: 3 p0 (3.75, 3.82, 3.90) x 1 seed per set, short chunks. Check: energy decreases, no stall-written
  `_EQ`, cell count, min cell area and signed-area orientation (collapse/inversion in QZ20), area histogram.
  Stop and report to the user if QZ20 cells collapse toward zero area.
- Full run with `batch_EQ_chunked.sbat` (CHUNK 1000, CAP 10000, USR1 save) -> residual-force check ->
  Hessian. Negative modes (double well): the LIBRARY (`calculateShearModulusWithEigenDecompositiond`) classifies each mode
  as positive (lambda > thr), zero (|lambda| <= thr, floppy test as now) or negative (lambda < -thr), and stores the
  lowest few eigenvalues and the negative-mode count in AGPerimTissue members. MainHess_claude.cpp only writes those values out;
  no threshold lives in a Main file. A negative mode is not ignored: push along the most negative mode and re-equilibrate, and
  record any tissue that stays stuck. Design: `vm_claude_plans/2026-10-21_negative-mode-push.md` (to discuss first).

## Step 4 - Bidisperse quadratic (2 ratios x 112 = 224 tissues)
- Inputs from the T2ss_OG raw inputs (`PhaseTransitionProject/T2ss_files/T2ss_OG/`; A0 = Abar = 1250, K 1.414214, Lambda 70.71;
  header and geometry identical to Aq2w, already QUADRATIC). Per tissue, choose 50
  big / 50 small cells at random (seeded from the tissue's rng), set per-cell A0 and P0 columns:
  r = 1.4: A0 = 1041.67 / 1458.33; r = 2: A0 = 833.33 / 1666.67; P0_c = p0*sqrt(A0_c). K = 1.4142 global, so each
  cell's energy is 1/2 K A0_c (a-1)^2 + 1/2 Lambda (P - P0_c)^2: every cell is a scaled copy (self-similar).
  Codes: BD14, BD20. Only the A0 and (new) P0 columns of the Cells section change.
- Pilot (3 p0 x 1 seed per ratio), then full run, residual forces, Hessian, G vs p0 against T2ss.
- **Bimodality check (pass/fail):** per tissue, histogram of A/Abar split by population; report
  mean(A_big)/mean(A_small) vs r and the overlap of the two populations (max small area vs min big area, and
  Ashman's D). Flag any tissue where they overlap.

## Bookkeeping
- Plans to `vm_claude_plans/2026-10-17_quartic2-and-bidisperse.md`; todo + lessons updated; README in the new
  work folder; back up to turbo (no `_error.txt`, no files > 100 MB); commit the library change and bookkeeping.

## Verification summary
Step 1 regression (uniform-A0 energies unchanged); Step 2 FD checks pass for all three C4 values and the
bidisperse tissue; pilots behave before full submission; residual forces of every new `_EQ` at the level of
the T2ss/QuOw references; bimodality report for bidisperse sets.
