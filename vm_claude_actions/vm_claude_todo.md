# Claude To-Do List

## PRIORITY 1
- [ ] ~~Write a standalone Python script...~~ **Dropped 2026-10-10 per user:** Python was only ever a quick side attempt; the shear modulus pipeline stays in C++ going forward (`AGPerimTissue.cpp`'s `calculateShearModulusWithEigenDecomposition(d)`), so no Python script is needed.

## Active Tasks
- [x] Fix Issue 3: unify eigenvalue thresholds across C++ single-prec and double-prec `calculateShearModulusWithEigenDecomposition(d)`; add floppy-state check (returns shear modulus = 0 for near-zero eigenvalue + non-negligible force projection). Also fixed `float currenteigenvalue` truncation bug in the double-precision function. Python script intentionally excluded — not used going forward (see lessons L002).
- [x] Fix Issue 5: unify normalization factor in both C++ shear-modulus functions — both now use `cells.size() / (K * pow(LX*LY, 2))`.
- [ ] ~~Fix Issue 4~~: dropped permanently, not deferred. `CalculateShearModulusForPythonCalculation` only existed to feed the Python path, which is no longer used — left untouched.
- [ ] Fix Issue 6: rename/verify `dA/dγ` functions (`GetdAdgammaForSimpleShearAlongXAxis` vs ALTERNATE); run finite-difference check to confirm which is correct. **Re-verified 2026-10-10: not a bug.** Docstrings in `AGVertexModel/AGCell.cpp` (~line 821 and ~936) explicitly mark the non-ALTERNATE version "THIS IS INCORRECT - IT DOESN'T HANDLE CORNER EDGES VERY WELL" and the ALTERNATE version "THIS IS THE CORRECT ONE!" — the code already calls the correct (ALTERNATE) one. Only `AGCell.h:107`'s comment is stale/confusing; no code change needed.
- [ ] Fix Issue 7: add finite-difference verification of `Getd2PdgammadvForSimpleShearAlongXAxis` against numerical Xi.
- [ ] Fix Issue 8: uncomment `d2Pdridrj_ALT - d2Pdridrj` debug print, run test, confirm agreement, then remove dead ALT computation.

## Backlog
- [ ] Review shear modulus and Hessian calculation code to understand current state (`AGStiffness`, `AGPerimModuli`, `AGPerimHessian_test`, `script_TO_CALC_SM_FROM_CPP.py`)

## Completed
- [x] Created initial `CLAUDE.md`
- [x] Created `claude_actions/` directory with `claude_todo.md`, `claude_recentactions.md`, `claude_lessons.md`
- [x] Updated `CLAUDE.md` with Boris Chernsky's workflow template + project-specific context
