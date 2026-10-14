# Plan: Exponential-energy validation of the Hessian shear modulus + G plots + ±ε benchmark test

**Status:** saved 2026-10-14, NOT yet approved — user will review next session (top of vm_claude_todo.md). Nothing executed.

## Context
The T2ss (quadratic energy) Hessian shear modulus is now proven correct: the code G equals an
independent FD reference to ~5e-9 after the three sign fixes. The remaining ~2% offset from the
squeeze benchmark is pure-vs-simple shear anisotropy plus a one-sided fit bias (README §5 of
the scratch folder). The user wants:
(1) a plot of G_simple (and the averaged G) added to their to-do to look at;
(2) **top priority:** check the Hessian G against the squeeze data for the non-quadratic energies.
    The code names come from `CompressionInducedFluidizationData/ExperimentalConditions_CodeNames.csv`:
    ExTs, L2Ex, M2Ex and similar sets;
(3) a test of the reflect-about-y=x idea (which gives the −ε branch).

Constraints (memory): no existing repo/data file is modified for testing. All work goes in
`/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/T2ss_SM_Hess_testing/` (or a sibling
folder). Runs go through SLURM, `--account=abmondal1`. At the end, `rsync -a --ignore-existing` to the turbo backup.

## Datasets (user-confirmed scope): Ex base sets + quadratic L2Ts, M2Ts and their LP variants
All have 7 seeds (rng 100–160), 16 p0 values (3.75–3.90), LX = LY = 353.553391 and 100 cells.
D = `/nfs/turbo/lsa-dkluben/abmondal/CompressionInducedFluidizationData`, V = `D/CIF_VertexModelSimulationInputOutpudEFiles`, H = `D/CIF_HomeFolderBackupCouldHaveDuplicates`.
| Code | energy | A0 | Lloyd | EQ inputs (candidate) | dE + benchmark |
|---|---|---|---|---|---|
| ExTs | EXP | 1250 | 10 | `H/ExTs_files/ExTs_EQ` | `H/ExTs_dE`, smdataExTs.mat |
| L2Ex | EXP | 625 | 10 | `V/L2Ex_files/L2Ex_EQ` | `V/L2Ex_dE`, smdataL2Ex.mat (its .m saves "…HP", so verify by refitting) |
| M2Ex | EXP | 2500 | 10 | `V/M2Ex_files/M2Ex_EQ` | `V/M2Ex_dE`, smdataM2Ex.mat |
| L2Ts | QUAD | 625 | 10 | `V/L2Ts_files/L2Ts_EQ` | `V/L2Ts_dE`, smdataL2Ts.mat |
| M2Ts | QUAD | 2500 | 10 | `V/M2Ts_files/M2Ts_EQ` | `V/M2Ts_dE`, smdataM2Ts.mat |
| L2Ts_LP | QUAD | 625 | 5 | `V/L2Ts_LP_files` or `H/L2Ts_LP_files_2.0` / `_redoCE` (ambiguous) | `V/L2Ts_LP_dE` (several smdata*.mat) |
| M2Ts_LP | QUAD | 2500 | 5 | `V/M2Ts_LP_files` or `D/M2Ts_LP_files` (ambiguous) | `V/M2Ts_LP_dE` (smdataM2TsHP.mat) |
**Input matching rule (handles the ambiguous LP folders):** the first line of every dE file is
the unsqueezed energy. A candidate EQ file is the right one only if its `finalEnergy` equals
that value to 6 digits. The benchmark G is always recomputed from the dE files with the MATLAB
fit (already validated against smdataT2ss.mat), and compared with whichever smdata*.mat
reproduces it. If no EQ candidate matches a dE set, report it rather than guess.
Copy the matched inputs into a scratch subfolder per code, with md5 provenance as for T2ss.

## Steps
1. **To-do bookkeeping** (`vm_claude_actions/vm_claude_todo.md`): new TOP item = the non-quadratic
   comparison. Add a USER item: "look at analysis/G_vs_p0.png (simple, pure, ½(simple+pure) vs squeeze)".
2. **Regenerate the T2ss G plot** with the fixed code (`results_fix123`) plus G_pure_H and
   ½(G_simple + G_pure_H) from `results_squeeze`. Extend `analysis/compare_to_dE.py`: add a
   results-dir argument, then merge `squeeze_per_seed.tsv`. Same MATLAB-style figure, G/K vs p0.
3. **Extend the FD reference to every pressure type** (`src/RefModel.h`):
   per-cell `ptype` field; the energy branch, `prefactors()` wA = dE/dA, and the energy formula
   mirror `AGPerimTissue::CalculateFreeEnergy` (E_area = ½ K a0 f(A/a0), with f = (x−1)², 2e^{−x},
   −2 ln(1+x), quartic, quartic2.0). Remove the "allQuadratic → reference invalid" guard in
   `MainHessFDCheck.cpp`, and store the type instead. The library's `GetPressure` and
   `GetNegativeDPressureDArea` (AGCell.cpp:39, :74) are thereby FD-checked as well.
4. **Run MainHessFDCheck + MainSqueezeCheck** (build_fixed, fixed snapshot) on all 7 sets,
   about 7 × 112 tissues, each run as a SLURM array.
   Generalize the sbat files to take a parameters-file argument (`FDparameters_<code>.txt`).
   Output to `results_<code>/` and `results_squeeze_<code>/`.
   Expected: G_code_d = G_FD (≤1e-8); benchmark energies reproduced to 6 digits; G_simple/G_squeeze
   ≈ 1.00–1.03 with the same anisotropy/fit factorization. Anything else means a real
   energy-type-specific error, to be localized with the existing term-by-term diagnostics.
5. **Generalize `compare_to_dE.py` / `compare_squeeze.py`** to take a code name (paths, smdata
   variable names `p0byrng<code>`, the file-name prefix). Make one plot per dataset plus a summary table.
6. **±ε test (reflection idea), no new simulations:** add the negative branch to
   `MainSqueezeCheck` (relax at fx = 1/(1+ε), fy = 1+ε, which is exactly what a y=x-reflected input
   would give under the existing code). Fit in s = ±ln(1+ε): G = (2c₂ − c₁)/(4·LX·LY). Report
   G_sym vs G_pure_H vs G_squeeze for T2ss and the Ex sets. If G_sym ≈ G_pure_H to ≲0.2%, the user can
   regenerate the real benchmark with reflected EQ files.
7. Update the scratch README (new sections), vm_claude_todo/recentactions, memory, and the rsync backup.

## Critical files (scratch unless noted)
- `src/RefModel.h`, `src/MainHessFDCheck.cpp`, `src/MainSqueezeCheck.cpp`, `Makefile`
- `batch_fdcheck.sbat`, `batch_squeeze.sbat` (parameters-file argument)
- `analysis/compare_to_dE.py`, `analysis/compare_squeeze.py`
- repo: `vm_claude_actions/vm_claude_todo.md`, `vm_claude_recentactions.md` (logs only, no source)

## Verification
- T2ss rerun after the RefModel change reproduces the previous summaries bit-for-bit (quadratic branch unchanged).
- Each Ex dataset: E_ref = CalculateFreeEnergy to ~1e-14; FD Hessian step stability ~1e-10;
  G_code_d vs G_FD; refit of the dE files reproduces smdata*.mat exactly.
- ±ε: G_sym from the frozen-topology energies vs G_pure_H, per seed.
