# Claude To-Do List

## PRIORITY 1
- [ ] Write a standalone Python script that computes the shear modulus from an equilibrated vertex model output file, implementing the full Hessian + affine-force-vector pipeline from scratch (no dependence on the C++ eigendecomposition output). Must follow the corrected LEBC formulas in `vertex_model_master_derivations.tex` and `vm_claude_proposed_fixes.md`.

## Active Tasks
- [ ] Fix Issue 3: unify eigenvalue thresholds across C++ single-prec, double-prec, and Python; add floppy-state check.
- [ ] Fix Issue 4: switch `CalculateShearModulusForPythonCalculation` to use `calculateHessiand` (double-precision Hessian export).
- [ ] Fix Issue 5: unify normalization factor in both C++ shear-modulus functions.
- [ ] Fix Issue 6: rename/verify `dA/dγ` functions (`GetdAdgammaForSimpleShearAlongXAxis` vs ALTERNATE); run finite-difference check to confirm which is correct.
- [ ] Fix Issue 7: add finite-difference verification of `Getd2PdgammadvForSimpleShearAlongXAxis` against numerical Xi.
- [ ] Fix Issue 8: uncomment `d2Pdridrj_ALT - d2Pdridrj` debug print, run test, confirm agreement, then remove dead ALT computation.

## Backlog
- [ ] Review shear modulus and Hessian calculation code to understand current state (`AGStiffness`, `AGPerimModuli`, `AGPerimHessian_test`, `script_TO_CALC_SM_FROM_CPP.py`)

## Completed
- [x] Created initial `CLAUDE.md`
- [x] Created `claude_actions/` directory with `claude_todo.md`, `claude_recentactions.md`, `claude_lessons.md`
- [x] Updated `CLAUDE.md` with Boris Chernsky's workflow template + project-specific context
