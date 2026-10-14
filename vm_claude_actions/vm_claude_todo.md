# Claude To-Do List

## PRIORITY 1
- [ ] ~~Write a standalone Python script...~~ **Dropped 2026-10-10 per user:** Python was only ever a quick side attempt; the shear modulus pipeline stays in C++ going forward (`AGPerimTissue.cpp`'s `calculateShearModulusWithEigenDecomposition(d)`), so no Python script is needed.

## Active Tasks
- [x] Fix Issue 3: unify eigenvalue thresholds across C++ single-prec and double-prec `calculateShearModulusWithEigenDecomposition(d)`; add floppy-state check (returns shear modulus = 0 for near-zero eigenvalue + non-negligible force projection). Also fixed `float currenteigenvalue` truncation bug in the double-precision function. Python script intentionally excluded — not used going forward (see lessons L002).
- [x] Fix Issue 5: unify normalization factor in both C++ shear-modulus functions — both now use `cells.size() / (K * pow(LX*LY, 2))`.
- [ ] ~~Fix Issue 4~~: dropped permanently, not deferred. `CalculateShearModulusForPythonCalculation` only existed to feed the Python path, which is no longer used — left untouched.
- [ ] Fix Issue 6: rename/verify `dA/dγ` functions (`GetdAdgammaForSimpleShearAlongXAxis` vs ALTERNATE); run finite-difference check to confirm which is correct. **Re-verified 2026-10-10: not a bug.** Docstrings in `AGVertexModel/AGCell.cpp` (~line 821 and ~936) explicitly mark the non-ALTERNATE version "THIS IS INCORRECT - IT DOESN'T HANDLE CORNER EDGES VERY WELL" and the ALTERNATE version "THIS IS THE CORRECT ONE!" — the code already calls the correct (ALTERNATE) one. Only `AGCell.h:107`'s comment is stale/confusing; no code change needed.
- [x] Fix Issue 7 — **diagnosed 2026-10-14** by the FD check (`/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/T2ss_SM_Hess_testing/README.md`). The perimeter parts of Ξ (incl. `Getd2PdgammadvForSimpleShearAlongXAxis`) are correct to 1e-10, and the sign convention is s = +1 throughout. The area part has bugs 2 and 3 below.
- [x] Fix Issue 8 — **diagnosed 2026-10-14**. The standard d2P/dridrj equals `d2Pdridrj_ALT` to 1e-12, and both match FD to 1e-11. ALT is redundant (dead-code removal optional, not done).
- [x] **Bug 1 (FD check 2026-10-14): Hessian d²A/dr dr block had the wrong sign** in `calculateHessian` and `calculateHessiand`. The "CCW convention" signs were applied to the CW vertex list. Fixed; `H_code` now matches the FD Hessian to 1.8e-11 (median, 128 tissues). The spurious soft modes and floppy returns at p0 3.80–3.82 are gone.
- [ ] **BOOKKEEPING FOR THE USER — correct the equations in `vertex_model_master_derivations.md`/`.tex`** (the code bug traced back to these):
  - [ ] §3.2 "Area Stress Term" ½[[0,1],[−1,0]] for j = i+1 is right **only when i → i+1 is counterclockwise**. Say so explicitly. The code stores vertices clockwise (`GetClockWiseVertices`), and applying this formula verbatim to the CW list caused bug 1. For a CW list the block is ½[[0,−1],[1,0]].
  - [ ] §2.1 the area gradient is internally inconsistent: ½ ẑ×(r_{i+1}−r_{i−1}) = ½(−Δy, Δx), but the explicit vector written is ½(Δy, −Δx). The explicit one is correct for CCW order.
  - [ ] §4.2 "∂A_α/∂γ = 0" is false per cell under LEBC. Boundary-crossing cells have dA/dγ ≠ 0 (only the sum over cells vanishes; the FD check confirms values ~10³). So §5 Ξ is missing the K_A ∇A_α ∂A_α/∂γ term, and §6.1 G_aff is missing K_A (∂A_α/∂γ)². The code includes both, correctly.
  - [ ] §4.1/§4.3/§5.1 use the homogeneous-affine-shear forms (all edges, x → x + γy), while §1.3/§8.4 use the Lagrangian LEBC form (only y-crossing edges shift by n_y γ L_y). The code uses LEBC. Pick one convention and write the strain derivatives consistently. In LEBC form: ∂P/∂γ = Σ_crossing n_y L_y Δx/l; ∂²A/∂γ∂y_j = σ (L_y/2)(ℓ_{j−1} − ℓ_{j+1}) with ℓ = cumulative image level; ∂²A/∂γ∂x_j = 0.
- [x] **Bug 2 (fixed 2026-10-14):** `AGCell::Getd2AdgammadyForSimpleShearAlongXAxis` had the wrong sign for corner cells (122/122 bad cells crossed both x and y). Rewritten from the shoelace formula on the unwrapped polygon: d²A/dγdy_j = σ (L_y/2)(ℓ_{j−1} − ℓ_{j+1}).
- [x] **Bug 3 (fixed 2026-10-14):** `AGCell::GetdAdgammaForSimpleShearAlongXAxisALTERNATE` flipped sign for 41/1546 y-crossing cells (the `crosscellvec.x` heuristic). Rewritten as dA/dγ = L_y Σ_k ℓ_k ∂A/∂x_k. This also revises Issue 6: ALTERNATE was not fully correct; it is now. The non-ALTERNATE version is unchanged and unused.
  - Verification on 128 tissues: 0 bad cells; Ξ matches FD to 1e-10. G_code_d = G_FD to ~5e-9 relative for every tissue with p0 ≤ 3.82. For p0 ≥ 3.83 the remaining differences are floppy-state returns in genuinely fluid tissues.
- [ ] **THEN:** investigate the median G_FD/G_squeeze = 1.04 in the solid (p0 ≤ 3.78); the user wants this closer to 1.

## Backlog
- [ ] Review shear modulus and Hessian calculation code to understand current state (`AGStiffness`, `AGPerimModuli`, `AGPerimHessian_test`, `script_TO_CALC_SM_FROM_CPP.py`)

## Completed
- [x] Created initial `CLAUDE.md`
- [x] Created `claude_actions/` directory with `claude_todo.md`, `claude_recentactions.md`, `claude_lessons.md`
- [x] Updated `CLAUDE.md` with Boris Chernsky's workflow template + project-specific context
