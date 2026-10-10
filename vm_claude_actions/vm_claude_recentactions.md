# Claude Recent Actions

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
