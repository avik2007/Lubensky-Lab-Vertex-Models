**SUPERSEDED 2026-10-26 (QUARTIC2.0 part):** the QZ05/QZ10/QZ20 one-A0-per-C4 design is replaced by A0 sweeps at fixed C4 (`2026-10-26_quartic2-A0-sweeps.md`). The bidisperse part (Step 4) and review #2 still stand.

**Status:** APPROVED 2026-10-22 (refreshed version; the 2026-10-17/21 draft is kept in git history). Two mandatory code reviews: #1 negative-mode probe (Step 0), #2 QUARTIC2.0 library fix (end of Step 1).

# Plan: QUARTIC2.0 validation + production, bidisperse (quadratic) tissues
(Refresh of `vm_claude_plans/2026-10-17_quartic2-and-bidisperse.md`, re-checked against the code on 2026-10-22.)

## Context
We need a regime where the quartic area term dominates, and separately want to know whether the model runs with a strongly
bidisperse tissue whose EQUILIBRIUM areas are bimodal. In the old QUARTIC form the quartic coefficient is K·C4, so K = 0 would
remove all area elasticity. The user's QUARTIC2.0 form
  E_A = ½ K A0 [(1−6C4)(a−1)² + C4(a−1)⁴],   a = A/A0,
has ½ f'' = 1 at |a−1| = 1 for every C4. It is pure quartic at C4 = 1/6 and a double well for C4 > 1/6.

### Code check (2026-10-22, `AG Vertex Model/AGVertexModel/`)
- **Forces and Hessian already support QUARTIC2.0 and per-cell A0.**
  - `AGCell::GetPressure` (AGVertexModel/AGCell.cpp:39) is the exact −dE/dA.
  - `AGCell::GetNegativeDPressureDArea` (:74) is the exact d²E/dA².
  - Both use the cell's own A0.
  - Note: AGCell.cpp has CRLF line endings and non-ASCII bytes, so plain `grep` skips it as binary; use `grep -a`.
- **Bug, confirmed:** `AGPerimTissue::CalculateFreeEnergy` (AGPerimVertex/AGPerimTissue.cpp:457) and `CalculatePressureEnergy` (:898)
  multiply the area sum by the LAST cell's `a0` after the loop (`0.5*K*a0*area_difference_sqr`). This is exact for uniform A0
  and wrong for bidisperse tissues. This energy drives the backtracking step test and the equilibrium slope test, so a mismatch
  with the forces can stall the step size and get written out as `_EQ`. `AGStiffTissue::CalculateFreeEnergy` does the per-cell
  sum correctly.
- **Gap, confirmed:** the Cells reader (:318–405) reads a per-cell `A0` column but sets `(*cptr).P0 = P0` (the header value) for
  every cell. Raw T2ss_OG files have the columns `id A0 EdgeNumb EdgeList Potential TCJStrength ThetaAn AnStrength
  CellPressureType` (value `Quadratic`, which falls through to the quadratic branch).
- **Done since the draft:** the negative-mode probe (`MainNegMode_claude.cpp`, plan `2026-10-22_negative-mode-probe.md`). It is
  implemented, awaiting your review, and replaces the draft's "store eigenvalues in members" idea. All QZ Hessian runs use it.

### Physics note for QZ20 (A0 = 2 Abar, C4 = 1/2), to watch in the pilot
- f(x) = (1−6C4)x² + C4x⁴ with x = a − 1 has minima at x² = (6C4−1)/(2C4) = 2, so a = 1 ± √2, i.e. **a = 2.41 or a = −0.41**.
- The left well is at a NEGATIVE area, and at the mean area a = 0.5 the area energy is concave (f'' ∝ −1.25).
- So the expected instability is: cells trade area, some grow toward a ≈ 2.4, and others shrink toward zero area or invert.
  `GetArea` is `fabs`, so an inverted cell reads as positive. With apoptosis OFF nothing removes collapsing cells.
- The QZ20 pilot gate below checks exactly this.

## User decisions (carried over)
- QUARTIC2.0 cases (K = 1.414214, the T2ss value, for all three; not K ∝ A0):
  | code | source raw set | A0/Abar | C4 | APOPTOSISON |
  |---|---|---|---|---|
  | QZ05 | Lq2w | 0.5 | 1/12 | 1 |
  | QZ10 | Aq2w | 1 | 1/6 | 1 |
  | QZ20 | Qu2w | 2 | 1/2 | 0 |
  - QZ10 sits at a ≈ 1, where pure quartic has ZERO area stiffness, so expect soft area modes.
  - K/A0 = 2.26e-3 / 1.13e-3 / 5.66e-4 differs across the sets.
- Bidisperse: quadratic only, from T2ss_OG raw files (header and geometry identical to Aq2w). 50/50 mix, mean A0 = Abar,
  big:small = 1.4:1 (BD14) and 2:1 (BD20), same shape index P0_c = p0·√A0_c. The equilibrated areas must come out
  near-bimodal (pass/fail).
- Step 1 mimics `AGStiffTissue::CalculateFreeEnergy`. Thresholds stay in the library. MANDATORY user code review after Step 1,
  before any pilot/production/commit of library code.

## Step 0 — CODE REVIEW #1: negative-mode probe (before it is used for anything, incl. QZ runs)
Hand the user the diff, then STOP until they approve:
- `AGPerimTissue.h/.cpp`: `REL_EIGENVALUE_THRESHOLD`, `FindNegativeModesd`, `PushAlongModed`, `PushAlongMostNegativeModed`
- `MainNegMode_claude.cpp`
- `NegMode/batch_negmode.sbat` and its Makefile
Present the smoke-test results (job 64634314) with it:
- **Regression:** Smd is bit-identical (probe = MainHess_claude = stored relthr) on all 11 T2ss/QuOw tissues.
- **Push mechanics on stiff modes:** at s = 1e-3, the curvature along a mode equals λ to 1e-5 to 1e-8 relative, and undoing a push
  restores the vertices to 1e-13. Gaps:
  - At s = 1e-4, round-off in E dominates, as expected.
  - 0 boundary crossings happened, so wrapping is still untested.
- **Detection:** 4 of 18 C4 ≠ 1 tissues have 1 negative mode, all fluid (p0 3.86). Also QuOw rng100 3.88 (2 modes) and T2ss
  rng150 3.88 (7 modes). Each has λ_min ≈ −1e-5, just past the threshold of ≈ 7e-6, inside a cluster of ~100 zero modes.
- **Flaw found:** for these weak modes the push cannot pick a sign.
  - ½|λ|s² ≈ 6e-12 is below one unit in the last place (ULP) of the total energy (E ≈ 3e5, ULP ≈ 6e-11).
  - It is also below the linear term from the residual force (|g·u|·s ~ 1e-7).
  - The soft-mode curvature check is off by 100×, for the same round-off reason.
  - So the outcomes NOPUSH ×3 and PUSHED ×1 are round-off noise. Round 2's doubled push (ΔE ≈ 2.4e-11) does not fix this.
- Decide in the review:
  - (a) whether λ ≈ −1e-5 in a fluid tissue is a real saddle or residual-force/round-off smear of the zero modes (i.e. whether the
    negative cut should be wider than the zero cut)
  - (b) how to choose the push sign when ΔE is below the energy resolution (e.g. per-cell energy differences, or the force
    projection)
  - (c) a test that exercises boundary wrapping
- No change to the probe until the user decides.

## Step 1 — Library fix (AGPerimTissue.cpp)
1. `CalculateFreeEnergy` and `CalculatePressureEnergy`: accumulate `a0 * f(a)` INSIDE the loop for every pressure type, and drop
   the trailing `a0` factor.
2. Cells reader: if the Cells section has a `P0` column, set `(*cptr).P0` from it; otherwise use the header P0 (old behavior).
   Check `HelperWriteOutputData_Write2File` (:1157): if it already writes a per-cell P0 column, restarts round-trip; if not, add
   one so `_EQ/_NE` files of bidisperse tissues keep P0_c.
3. Verbose comments; add a section to `vm_claude_reports/2026-10-15_vertex_model_code_changes_report.md`.
- Regression: on uniform-A0 files (T2ss, QuOw, new QuTs/MTQ3 EQ) the new energy equals the stored `finalEnergy` to ~1e-14
  relative. Probe Smd is unchanged (bit-identical expected, since the Hessian does not use these functions).
- Build in a NEW work folder `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Quartic2_Bidisperse/` (own binaries:
  equilibration_chunk, hess, negmode). Never rebuild the binaries the running Quartic_C4fix jobs use (L007).
- **CODE REVIEW #2: QUARTIC2.0 / per-cell energy + P0 reader.** Hand the user the Step 1 diff plus the regression numbers, then
  STOP. No Step 2 test runs, pilots, production or commit of library code until the user approves.

## Step 2 — Verification (before any production)
- Extend the independent checks to QUARTIC2.0 and per-cell A0/P0:
  - `Quartic_C4fix/analysis/residual_forces/vmcheck.py` (energy, forces)
  - `Quartic_C4fix/src/RefModel.h` + `MainHessFDCheck.cpp` (Hessian, Ξ, G vs finite differences)
- Test tissues: Lq2w/Aq2w/Qu2w raw and EQ states converted to QUARTIC2.0 at C4 = 1/12, 1/6, 1/2, plus one bidisperse tissue.
- Pass criteria:
  - library energy = reference
  - force = −FD gradient (1e-9 relative)
  - Hessian = FD Hessian (~1e-10)
  - G_code = G_FD on non-floppy states
  - For QZ20 also check that the probe's λ_min and its push behave (expected to show real negative modes; this is the
    end-to-end test the negative-mode smoke test could not do).

## Step 3 — QUARTIC2.0 production (3 × 112 = 336 tissues)
- Inputs: copy the source raw sets (same geometries, rng 100–160, p0 3.75–3.90). A Python script (LF endings, exact string edits)
  changes only:
  - the header C4 (full precision), K = 1.414214, APOPTOSISON
  - the cell type `QUARTIC` → `QUARTIC2.0`
- Add the new codes as rows in `PhaseTransitionProject/ExperimentalConditions_CodeNames.csv` (no name clashes).
- Pilot: 3 p0 (3.75, 3.82, 3.90) × 1 seed per set, chunked. Check:
  - the energy decreases
  - no `_EQ` written because the step stalled
  - cell count
  - **min cell area and signed-area orientation** (QZ20 collapse/inversion, see the physics note)
  - area histogram
  - **Stop and report if QZ20 cells collapse toward zero area or invert.**
- Full run: `batch_EQ_chunked.sbat` (CHUNK 1000, CAP 20000, USR1 save) → residual-force check → Hessian via
  `batch_negmode.sbat` / `MainNegMode_claude` (negative-mode ledger per tissue; STUCK handling decided from this data).

## Step 4 — Bidisperse quadratic (2 × 112 = 224 tissues)
- From T2ss_OG raw inputs (A0 = Abar = 1250, K 1.414214, Λ 70.71). Per tissue, choose 50 big / 50 small cells at random, seeded
  from the tissue's rng.
  - r = 1.4: A0 = 1041.67 / 1458.33. r = 2: A0 = 833.33 / 1666.67.
  - Add a per-cell P0 column with P0_c = p0·√A0_c. Only the A0 and P0 columns change.
- Pilot (3 p0 × 1 seed per ratio) → full run → residual forces → Hessian (negmode driver) → G vs p0 against T2ss.
- **Bimodality check (pass/fail):** per tissue, A/Abar split by population; mean(A_big)/mean(A_small) vs r; overlap (max small vs
  min big, Ashman's D). Flag any tissue where the populations overlap.

## Verification summary
- Step 1: uniform-A0 energies unchanged (regression).
- Step 2: FD checks pass for all three C4 values and the bidisperse tissue.
- Pilots behave before the full submission.
- Residual forces of every new `_EQ` at the level of the T2ss/QuOw references.
- Bimodality report for the bidisperse sets.
- Bookkeeping: plan file in `vm_claude_plans/` (update the 10-17 file and mark it approved), todo, README in the work folder,
  turbo backup (no `_error.txt`, no files > 100 MB), commit after the user's review.
