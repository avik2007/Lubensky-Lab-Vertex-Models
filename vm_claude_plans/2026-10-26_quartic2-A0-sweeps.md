**Status:** APPROVED by the user 2026-10-26 (verbally: "submit the quartic 2.0"). C4 = 1/6 sweep submitted (job 64834279); the 18-tissue pilot was skipped at the user's request and replaced by monitoring the first chunks.

# Plan: QUARTIC2.0 A0 sweeps — C4 = 1/12 first, then C4 = 1/6 (QZ20 later)

## Context
The approved QUARTIC2.0 plan (`vm_claude_plans/2026-10-17_quartic2-and-bidisperse.md`) paired each C4 with ONE A0:
QZ05 (0.5Ā, C4 = 1/12), QZ10 (Ā, 1/6), QZ20 (2Ā, 1/2). None of it has run. The user wants it rebuilt as A0 sweeps at fixed C4,
like every other figure (L2Ts/T2ss/M2Ts; LqOw/AqOw/QuOw; Lq2w/Aq2w/Qu2w). This is the core data for the paper.
User decisions (2026-10-26):
- Order: C4 = 1/6 (pure quartic) sweep first, then C4 = 1/12 (user changed the order 2026-10-26), then QZ20 (double well; still needs the negative-mode probe review #1).
- A0 = 0.5Ā, Ā, 2Ā; **K ∝ A0** (0.707 / 1.414 / 2.83, the triplet convention, so K/A0 is the same in all three).
  (Explained to the user: K scales both area terms, so it does not change the quartic:quadratic ratio C4·x²/(1−6C4).)
- High-polydispersity (_LP) versions: later, after these two sweeps.
- Capped tissues (NE at CAP) are left out of the G analysis; dropped counts reported per (set, p0).

Energy: E_A = ½·K·A0·[(1−6C4)x² + C4x⁴], x = A/A0 − 1. Both C4 values are convex (C4 ≤ 1/6). At C4 = 1/6 and A0 = Ā the area
stiffness at x = 0 is zero, so soft area modes and slow equilibration are expected there.

## Sets (6 × 112 = 672 tissues; rng 100–160, p0 3.75–3.90)
| code (proposed) | C4 | A0 | K | source raw set (same starting tissues for all) |
|---|---|---|---|---|
| L2Z1 | 1/12 | 0.5Ā = 625 | 0.707107 | Lq2w raw |
| AqZ1 | 1/12 | Ā = 1250 | 1.414214 | Aq2w raw |
| M2Z1 | 1/12 | 2Ā = 2500 | 2.828427 | Qu2w raw |
| L2Z6, AqZ6, M2Z6 | 1/6 | same three | same | same |
The source headers already carry these A0 and K values. Only two things change: the header C4 (written at full precision) and the
cell type `QUARTIC` → `QUARTIC2.0`. APOPTOSISON stays 1. Check the codes against `PhaseTransitionProject/ExperimentalConditions_CodeNames.csv`
before use (no clash today with QZ*), then add rows.

## Code: committed library only (no review gate)
- Uniform A0 per tissue, so the HEAD energy is exact (the per-cell a0 fix in review #2 matters only for bidisperse tissues).
  HEAD `CalculateFreeEnergy` already has the QUARTIC2.0 branch. `AGCell::GetPressure` and `GetNegativeDPressureDArea` are exact.
- New work folder `claude_VM_runs/Quartic2_Sweep/` with its OWN binaries, built from HEAD library sources
  (`git stash`-free: build from a `git worktree` of HEAD, so the uncommitted review-#1/#2 edits are not compiled in):
  `agperim_equilibration_claude_chunk` (MainAGPerim_claude.cpp), `agperim_hess_claude` (MainHess_claude.cpp). Copy
  `Quartic_C4fix/batch_EQ_chunked.sbat`, `Makefile` and `batch_SM_relthr.sbat`. Never touch the Quartic_C4fix binaries
  (jobs 64823317 and 64745543 are using them; lesson L007).

## Steps
1. **Inputs**: Python script (`-I`, exact string edits, keep CRLF/LF as found) writes the 336 C4 = 1/12 files plus the 336 C4 = 1/6 files to
   `Quartic2_Sweep/raw/<code>/`. Check that each file differs from its source only in the C4 line and the cell-type column (diff count).
2. **FD check (Step 2 of the old plan, reduced to uniform A0)**: extend the independent reference
   (`Quartic_C4fix/src/RefModel.h` energy/gradient/curvature lines 47–203; `analysis/residual_forces/vmcheck.py` lines 8, 81–83) with a
   QUARTIC2.0 branch: E = ½KA0[(1−6C4)x² + C4x⁴], dE/dA = K[(1−6C4)x + 2C4x³], d²E/dA² = (K/A0)[(1−6C4) + 6C4x²].
   Run `MainHessFDCheck.cpp` on one raw and one short-evolved tissue per set (6 sets).
   Pass: energy = reference; force = −FD gradient to 1e-9 rel; Hessian = FD to ~1e-10; G_code = G_FD on non-floppy states.
3. **Pilot**: p0 = 3.75, 3.82, 3.90 × seed 100 × 6 sets = 18 tissues, chunked (CHUNK 1000, CAP 20000).
   Check: energy decreases; cell count stays 100; min area > 0; area histogram (A/A0 per set); EQ reached or a clean slope.
   Report to the user before the full submission.
4. **Production, C4 = 1/12 first**: `sbatch --array=1-336 --export=ALL,PARAMS=EQparameters_Z1.txt,CAP=20000 batch_EQ_chunked.sbat`.
   Resubmit SIGNAL stops with the same command. Submit C4 = 1/6 (`EQparameters_Z6.txt`) right after, if the user agrees (it is
   independent; the cluster can run both).
5. **Analysis per sweep**: Hessian G on every `_EQ.txt` (`agperim_hess_claude` → `relthr/<code>/`). For capped tissues, run the force/Hessian check
   (`ForceHessCheck/batch_forcehess.sbat`) and leave them out of G. Then one triplet figure per C4 in the Hess_SM_groups style
   (0.5Ā black o, 2Ā blue s, Ā red D), overlaid on L2Ts/T2ss/M2Ts. The G table: mean ± SEM (n) per p0, plus dropped counts.
6. **Then**: _LP versions of both sweeps; QZ20 (needs review #1 / negative-mode probe).

## Verification
- Input diff: only C4 and the type column change (all 672 files).
- FD check passes on all 6 sets before the pilot.
- Pilot report (energy, cell count, min area, areas) before production.
- G regression: `agperim_hess_claude` built from HEAD gives Smd identical to the stored relthr value on 3 T2ss tissues.
- Bookkeeping: supersede the QZ table in the 10-17 plan (note it there), todo, recent actions, `vm_file_structure.md` if repo files
  change (RefModel/vmcheck are in scratch, not the repo), README in `Quartic2_Sweep/`, turbo backup at session end.
