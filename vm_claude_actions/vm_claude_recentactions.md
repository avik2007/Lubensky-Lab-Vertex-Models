# Claude Recent Actions

## 2026-10-15 (Hessian G for the low-poly triplets)
- Added `AG Vertex Model/AGVertexModel/AGPerimHessian_test/AGPerimHessian_test/MainHess_claude.cpp`: copy of MainHessianTest.cpp with minimal changes (argv input, no evolution per user, debug dumps and the Python-path call removed, Smd now uses the double-precision function, per-cell CW vertex lists set, writes `<tissue>_SM.txt`). Listed in `vm_file_structure.md`.
- Built in `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Hess_SM_groups/` against the repo library in place; smoke test job 64086589 reproduced the FD-verified G_code_d to 5e-13 (the first smoke run, without the CW-vertex set-up, was ~96x off).
- SLURM jobs 64086617/8/9: 1007/1007 tissues completed. Inputs from `Avik/PhaseTransitionProject` (read-only). Plots `analysis/G_vs_shapeindex_{Ts,Lo,Ex}_group.png`; tables `G_hess_all.tsv`, `G_groups_summary.tsv`.
- Result: Hessian/squeeze ratio of seed means for p0 <= 3.80 is 1.01-1.04 in every dataset (L2Ex 0.84 at p0 = 3.80 only); both ~0 beyond 3.81.

## 2026-10-14 (squeeze-benchmark offset)
- Wrote `src/MainSqueezeCheck.cpp` and `analysis/compare_squeeze.py` in the scratch FD-check folder (`/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/T2ss_SM_Hess_testing/`); no repo source touched. SLURM job 63973299.
- Result: the Hessian vs squeeze-fit offset is sample anisotropy (pure vs simple shear) plus a one-sided finite-strain fit bias. It is not a code error. Details are in that folder's README §5; the next step is a user decision (see todo).
- Saved the next plan, unapproved, to `vm_claude_plans/2026-10-14_nonquadratic-and-LP-hessian-vs-squeeze-comparison-plus-pm-eps-test.md`. It is PRIORITY 0 in the todo.
- Close-out: synced the scratch folder to `/nfs/turbo/lsa-dkluben/abmondal/CompressionInducedFluidizationData/claude_hessian_testing/T2ss_SM_Hess_testing/` (checksum-verified). The 4 files that changed since the first backup (README.md, Makefile, batch_fdcheck.sbat, src/MainHessFDCheck.cpp) were updated there; their earlier versions are kept as `*.pre-final-2026-10-14`.

## 2026-10-10
- Picked up the approved-but-unimplemented plan from `vm_claude_plans/2026-10-09_shear-modulus-eigenvalue-threshold-normfactor-and-hessian-precision-fixes.md`.
- Independently re-verified its claims against current source before trusting them, including the one that looked most likely to be wrong (the `GetdAdgammaForSimpleShearAlongXAxis` vs `...ALTERNATE` naming claim) — confirmed accurate via `grep -a` (the file's non-standard encoding silently breaks plain `grep`, which initially made it look like neither function was implemented anywhere; `-a` revealed both, with docstrings settling which one is correct).
- User corrected scope: "ignore Python" means ignore anything that exists only to serve the Python path, not just the `.py` file — so Fix 3 from the 2026-10-09 plan (double-precision Hessian export for `CalculateShearModulusForPythonCalculation`) was dropped entirely rather than implemented. Recorded as lesson L002.
- Implemented in `AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp`:
  - **Issue 3**: both `calculateShearModulusWithEigenDecomposition` and `calculateShearModulusWithEigenDecompositiond` now share an `EIGENVALUE_THRESHOLD = 1e-5` / `FLOPPY_TOLERANCE = 1e-4` pattern — modes below threshold with non-negligible force projection (`xi_dot_u`) now return shear modulus = 0 (floppy-state early-return) instead of silently contributing (single-precision function) or being silently skipped (double-precision function). Also fixed `calculateShearModulusWithEigenDecompositiond`'s `currenteigenvalue` from `float` to `double` (was truncating before the threshold comparison).
  - **Issue 5**: `calculateShearModulusWithEigenDecomposition`'s normalization factor changed from `1/(K*A0^2*Nc)` to `Nc/(K*(LX*LY)^2)`, matching the double-precision function.
- Verified by compiling (`g++ -std=c++14 -fopenmp -fno-fast-math -Wall`) `MainHessianTest.cpp` + all its dependencies against the NuGet-vendored Eigen at `AGPerimHessian_test/packages/Eigen3.3.3.9/...` — full build and link succeeded (exit 0), only pre-existing warnings (unrelated unused variables etc.), no new errors.
- Updated `vm_claude_actions/vm_claude_todo.md`: marked Issues 3 and 5 done, Issue 4 dropped (not deferred), Issue 6 re-confirmed not-a-bug, dropped the standalone-Python-script Priority-1 item.
- **Committed** on `MondalShearDev`: commit `5f7a19504fcca4683762497d01f54911c0aef881` ("Fix shear modulus eigenvalue threshold, floppy-state, and normalization bugs"), 2026-10-10 21:54:46 -0400. Scope: exactly 4 files — `AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp` and the three `vm_claude_actions/*.md` files. Staged explicitly by path (not `git add -A`) because the working tree has ~2900+ unrelated dirty entries from pre-existing CRLF line-ending drift across the repo (predates this session, not touched). To revert just this fix: `git revert 5f7a19504fcca4683762497d01f54911c0aef881`. To inspect it later: `git show 5f7a19504fcca4683762497d01f54911c0aef881`.

## 2026-03-26 (session 2)
- Read `vertex_model_master_derivations.tex` in full.
- Explored all relevant C++ source files: `AGPerimTissue.cpp` (Hessian + shear modulus), `AGPerimCell.cpp` (LEBC strain derivatives), `AGCell.cpp` (area gradient, pressure), `EdgeClass.cpp` (`GetVectorDifference` LEBC logic).
- Compared document formulas against code implementations and produced `vm_claude_proposed_fixes.md` at repo root documenting 8 discrepancies:
  1. `dP/dγ` formula in document uses vertex-shear convention, not LEBC (code is correct).
  2. `d²P/dγ²` formula in document uses `Δy` instead of `Lᵧ` (code is correct).
  3. Eigenvalue thresholds inconsistent: C++ single-prec has none, double-prec uses 1e-5, Python uses 1e-4; document specifies 1e-14 with floppy-state early-return not implemented anywhere.
  4. `CalculateShearModulusForPythonCalculation` exports float32 Hessian — precision loss near the rigidity transition.
  5. Normalization factor differs between single- and double-precision shear-modulus functions.
  6. `dA/dγ` function names inverted: ALTERNATE function is actually used; standard one is discarded.
  7. `d²P/(dγ·dr_i)` mixed derivative unverified by finite difference.
  8. Edge-based Hessian ALT (`d2Pdridrj_ALT`) computed but silently discarded; comparison was commented out.
- Note: the `.tex` file was updated by the user during this session to correct the strain-derivative section and the shear-modulus prefactor (now `1/(4Nc)`).

## 2026-03-26 (session 1)
- Created `vm_file_structure.md` at repository root: full annotated directory tree of the entire repo, including descriptions of every major folder and file, the C++ class hierarchy, data flow diagram, and environment setup info.
- Updated `CLAUDE.md`: removed "Project Overview" subsections (C++ class hierarchy, data flow, environment setup) that duplicated structure info, replaced with a pointer to `vm_file_structure.md` and an explicit instruction to keep that file updated whenever files/directories are added or removed.

## 2026-03-19
- Created `CLAUDE.md` at repository root with project overview, environment setup, C++ class hierarchy, MATLAB/Python tool descriptions, and data flow summary.
- Created `claude_actions/` directory with `claude_todo.md`, `claude_recentactions.md`, and `claude_lessons.md`.
- Rewrote `CLAUDE.md` using Boris Chernsky's workflow template as the core, adding: scientific workflow context (prompts come from external LLM), branch scope restriction to shear modulus/Hessian work only, task management conventions pointing to `claude_actions/` files.
- Recorded L001 in `claude_lessons.md`: don't assume you know a named reference — ask first.
