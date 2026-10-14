# Fix shear modulus bugs in AGPerimTissue.cpp (C++ only)

**Status:** Approved by user 2026-10-09. Not yet implemented — user logged off
before coding started. Pick this up as the next task on `MondalShearDev`.

## Context

A previous session audited the shear-modulus/Hessian code and left an unactioned
list of discrepancies in `vm_claude_proposed_fixes.md` (repo root). I re-verified
every claim directly against the current source (`AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp`
and `AGVertexModel/AGCell.cpp`) rather than trusting the writeup blindly, and found:

- Some of the old findings are confirmed real bugs.
- One (naming of `GetdAdgammaForSimpleShearAlongXAxis` vs `...ALTERNATE`) is **not**
  a bug — the code already calls the function that its own docstring in `AGCell.cpp`
  marks `"THIS IS THE CORRECT ONE!"`, while the unused one is explicitly marked
  `"THIS IS INCORRECT"`. Only the header comment (`AGCell.h:107`) is stale/confusing.
  No code change needed there.
- One new issue was found while reading `calculateShearModulusWithEigenDecompositiond`:
  it declares `currenteigenvalue` as `float` even though the whole function is
  double-precision (`MatrixXd`), silently truncating the eigenvalue before the
  threshold comparison.

The user confirmed scope: fix only the three C++ issues below (not the Python
script, which is a separate, not-currently-trusted path). All changes are confined
to `AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp`, which is inside
the branch's stated scope (`AGPerimVertex` holds the Hessian/shear-modulus logic
that `AGStiffness`, `AGPerimModuli`, and `AGPerimHessian_test` all link against).

## Fix 1 — Unify eigenvalue threshold + add floppy-state check

Affects both `calculateShearModulusWithEigenDecomposition` (lines ~2739-2950,
single-precision/`MatrixXf`) and `calculateShearModulusWithEigenDecompositiond`
(lines ~2951-3168, double-precision/`MatrixXd`).

Current state:
- Single-precision function (loop at ~2898-2922): accumulates
  `1/eigenvalue * (Xi·u)^2` into the nonaffine correction for **every** mode with
  no threshold — a near-zero eigenvalue blows up the result.
- Double-precision function (loop at ~3112-3139): already guards with
  `if (currenteigenvalue > 1e-5)`, but `currenteigenvalue` is declared `float`
  (line 3016), truncating the double-precision eigenvalue before the comparison.
- Neither function implements the "floppy state" case: a near-zero eigenvalue
  whose mode still carries a non-negligible force projection (`Xi·u`) signals a
  mechanically unstable/marginal configuration where the quasi-static nonaffine
  correction diverges. The correct behavior (per the master derivation doc and
  Manning/Merkel's method) is to detect this and return shear modulus = 0 rather
  than a spurious finite number.

Plan:
- In the double-precision function: change `float currenteigenvalue;` (line 3016)
  to `double currenteigenvalue;`.
- In both functions, restructure the per-mode loop to:
  1. Compute `xi_dot_u = (eigenvectors().col(index).transpose() * FullBoundaryVector)(0,0)`
     once (currently recomputed inline via `pow(...)`).
  2. If `currenteigenvalue > EIGENVALUE_THRESHOLD`: accumulate
     `nonaffine_shearportion += xi_dot_u*xi_dot_u / currenteigenvalue` (existing
     behavior, just reusing `xi_dot_u`). Keep all existing debug `*errorfile <<`
     logging as-is.
  3. Else (eigenvalue at/below threshold): if `fabs(xi_dot_u) > FLOPPY_TOLERANCE`,
     log a clear message to `errorfile` identifying the floppy mode (index,
     eigenvalue, `xi_dot_u`) and `return 0.0` immediately from the function.
     Otherwise (genuinely negligible projection), skip the mode silently — same
     as today's double-precision behavior.
- Use the same two named constants in both functions:
  `EIGENVALUE_THRESHOLD = 1e-5` (reuses the value already proven in the existing
  double-precision code, rather than inventing an untested number — safe for both
  the float32 and float64 code paths) and `FLOPPY_TOLERANCE = 1e-4` (chosen with
  margin above float32 rounding noise). Document the reasoning in a comment at
  each definition site, and cross-reference `vm_claude_proposed_fixes.md` Issue 3.

## Fix 2 — Normalization factor in the single-precision function only

Verified: only `calculateShearModulusWithEigenDecomposition` (line 2927) disagrees.
`calculateShearModulusWithEigenDecompositiond` (line 3145) and
`CalculateShearModulusForPythonCalculation` (line 3356) already both use
`cells.size() / (K * pow(LX*LY, 2))`.

Plan: in `calculateShearModulusWithEigenDecomposition`, replace

```cpp
double normfactor =  1 / (K * pow((*this).A0, 2) * (*this).cells.size());
```

with

```cpp
// Canonical normalization — matches calculateShearModulusWithEigenDecompositiond()
// and CalculateShearModulusForPythonCalculation(): g* = Nc / (K * A_box^2), with
// A_box = LX*LY the actual simulation box area. This holds regardless of whether
// the tissue's target area (*this).A0 matches the true mean per-cell area.
double normfactor = (*this).cells.size() / (K * pow(LX * LY, 2));

// --- Old (pre-fix) normalization, kept for reference ---
// double normfactor_OLD = 1.0 / (K * pow((*this).A0, 2) * (*this).cells.size());
// Implicitly assumed (*this).A0 == (LX*LY) / Nc, i.e. the tissue's target area
// equals the actual average per-cell box area. That assumption breaks once A0
// is a fixed target that no longer matches the true average cell area (e.g.
// after cell division or anisotropic growth), producing a shear modulus that
// disagreed with the double-precision and Python-export functions for the same
// configuration. See vm_claude_proposed_fixes.md, Issue 5.
```

## Fix 3 — Double-precision Hessian export for Python

`CalculateShearModulusForPythonCalculation` (lines 3169-3376) exists solely to
hand off `Hessian`, `FullBoundaryVector`, `normfactor`, and `d2Edgamma2` to the
Python script as text (its own eigendecomposition loop is entirely commented
out — it does not compute a shear modulus itself). Right now every Eigen object
in it is `MatrixXf` (float32), including the Hessian and boundary vector that
Python will eigendecompose, defeating the purpose of a "more accurate Python
path."

Plan: in this function only, do a mechanical type swap:
- `calculateHessian(pythonfile)` → `calculateHessiand(pythonfile)` (the existing
  double-precision Hessian builder, already used elsewhere in the file).
- Change every local `MatrixXf` declaration to `MatrixXd`: `Hessian`,
  `hessianeigen` (`SelfAdjointEigenSolver<MatrixXd>`), `lu_decomp`
  (`FullPivLU<MatrixXd>`), `FullBoundaryVector`, `d2Edgamma2_mat`,
  `d2AdrdgammaVEC`, `d2PdrdgammaVEC`, `dAdrdAdgammaVEC`, `dPdrdPdgammaVEC`, and the
  per-cell `BoundaryVector` declared inside the cell loop (~line 3241).
- Change `float currenteigenvalue;` (line 3230) to `double` for consistency
  (it's currently unused since that loop is commented out, but the function
  should be double-precision throughout).
- No change to output formatting, markers, or what gets written to `pythonfile` —
  the Python script's text parsing (`script_TO_CALC_SM_FROM_CPP.py`) keys off
  fixed string markers (`"normfactor = \n"`, `"Hessian = \n"`, etc.) that are
  untouched by this type change, just populated with higher-precision digits.

## Out of scope (confirmed, not touched)

- `script_TO_CALC_SM_FROM_CPP.py` — user explicitly deferred this. Separately
  noted during this audit (not in `vm_claude_proposed_fixes.md`): it calls
  `np.linalg.eig` instead of `np.linalg.eigh` on a symmetric Hessian, which can
  return complex-dtype eigenvalues/eigenvectors and silently fail inside the
  script's bare `except:` block. Worth fixing later but not requested now.
- `GetdAdgammaForSimpleShearAlongXAxis` vs `...ALTERNATE` naming — confirmed not
  a bug; leaving code as-is. (Could optionally fix the stale header comment in
  `AGCell.h:107` in a future pass — not requested now.)
- Issues 7 (mixed derivative sign convention) and 8 (unused `d2Pdridrj_ALT`) from
  `vm_claude_proposed_fixes.md` — these need finite-difference verification before
  any code change, which is a separate task (matches the existing Priority-1 item
  in `vm_claude_todo.md`), not requested in this round.

## Verification

- This is a Visual Studio C++ project (`.vcxproj`/`.sln`) targeting Eigen 3.3.9 via
  NuGet, built on Windows normally. Check whether a Linux-side build path exists
  (the repo has "Great Lakes/local makefiles" per the latest commit) and attempt to
  compile `AGPerimVertex`/`AGStiffness` with it to catch type errors from the
  `MatrixXf`→`MatrixXd` swap.
- If no build is feasible in this environment, do a careful manual re-read of the
  diff to confirm every `MatrixXf` touched in Fix 3 was changed consistently (no
  leftover float/double mismatch that would fail to compile), and that both
  threshold loops in Fix 1 compile-check by inspection (matching brace structure,
  `xi_dot_u` declared once per loop iteration, no shadowing).
- Update `vm_claude_actions/vm_claude_todo.md` and `vm_claude_actions/vm_claude_recentactions.md`
  to mark Issues 3, 4, 5 as fixed and record what was done, per this repo's own
  CLAUDE.md task-tracking convention.
