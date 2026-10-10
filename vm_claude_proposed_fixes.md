# Proposed Fixes: Hessian and Shear Modulus Implementation

This document identifies discrepancies between the canonical derivation document
(`vertex_model_master_derivations.tex`) and the C++ / Python implementation.
For each issue the document formula, the code formula, a physical analysis, and a
proposed resolution are given.

All C++ line numbers refer to
`AG Vertex Model/AGVertexModel/AGPerimVertex/AGPerimTissue.cpp` unless noted otherwise.

---

## Issue 1 — `dP/dγ` formula: document uses vertex-shear convention, not LEBC

### What the document says (§4)

```
dP_α/dγ = Σ_{<ij> ∈ α}  Δx_ij · Δy_ij / l_ij
```

where the sum is over **all edges** of cell α, and `Δx_ij`, `Δy_ij` are the
(unwrapped) edge displacements.

### What the code does

`AGPerimCell::GetdPdgammaForSimpleShearAlongXAxis()` (~line 350, `AGPerimCell.cpp`):

```cpp
// Only y-crossing edges (yflag == 1) contribute:
directioncos = edgeVect.x / edgeVect.Magnitude();   // Δx/l
updownfactor = (edgeVect.y < 0) ? -1 : 1;           // sign of Δy_PBC
dPdgamma += LY * updownfactor * directioncos;        // LY · sign(Δy_PBC) · Δx/l
```

The sum runs over **y-crossing edges only** and uses `LY` instead of `Δy_PBC`.

### Why they differ — which is correct

Under LEBC, changing the shear parameter γ by δγ shifts the top periodic image by
`δγ · LY` in x. For a y-crossing edge with winding number `n_y = sign(Δy_PBC)`:

```
d l / dγ = n_y · LY · Δx_phys / l
```

This is exactly what the code computes. The document formula `Δx · Δy / l` is the
"vertex shear" formula valid when every vertex position changes affinely with γ. Under
LEBC only the boundary condition shifts; vertex positions are free to relax, so the
vertex-shear formula gives different numerical values (off by the factor `LY / |Δy_PBC|`
per crossing edge).

**The code is physically correct for LEBC. The document formula is incorrect for LEBC.**

### Proposed fix

Update §4 of `vertex_model_master_derivations.tex` to the LEBC formula:

```
dP_α/dγ = Σ_{<ij> crossing y-boundary} n_ij · L_y · Δx_ij / l_ij
```

where `n_ij = sign(Δy_ij^PBC)` is the y-winding number of the edge
(+1 for edges that exit through the top boundary, −1 for the bottom).

---

## Issue 2 — `d²P/dγ²` formula: document and code give different expressions

### What the document says (§4)

```
d²P_α/dγ² = Σ_{<ij> ∈ α}  Δy²_ij / l_ij · (1 − Δx²_ij / l²_ij)
```

Expanding `(1 − Δx²/l²) = Δy²/l²`, this becomes `Σ Δy⁴/l³` when `Δy` is
interpreted as the PBC minimum-image y-displacement.

### What the code does

`AGPerimCell::Getd2Pdgamma2ForSimpleShearAlongXAxis()` (~line 399, `AGPerimCell.cpp`):

```cpp
directionsin = abs(evec.y) / edgelength;          // |Δy_PBC| / l
d2Pdgamma2  += LY * LY * pow(directionsin, 2) / edgelength;
// = LY² · Δy²_PBC / l³
```

### Why they differ — which is correct

Taking the second derivative of `l(γ) = sqrt((Δx + n_y · γ · LY)² + Δy²)`:

```
d²l/dγ² = LY² / l  −  (n_y · LY · Δx)² / l³
         = LY² / l · (1 − Δx² / l²)
         = LY² · Δy²_PBC / l³
```

This is what the code computes and it is correct for LEBC. The document expression
`Δy²/l · (1−Δx²/l²)` equals the code expression only when `|Δy_PBC| = LY`, i.e.,
all crossing-edge vertices lie exactly on y = 0 or y = L_y — not the general case.

**The code is physically correct. The document formula replaces `LY` with `Δy_PBC`
in error.**

### Proposed fix

Update §4 of `vertex_model_master_derivations.tex`:

```
d²P_α/dγ² = Σ_{<ij> crossing y-boundary}  L_y² · Δy²_ij / l³_ij
           = Σ  L_y² / l_ij · (1 − Δx²_ij / l²_ij)
```

---

## Issue 3 — Eigenvalue threshold inconsistency across all implementations

### What the document says (§7 Algorithm)

```
If ω²_m > 10⁻¹⁴:  add (Ξ · u_m)² / ω²_m to the nonaffine correction
If ω²_m < 10⁻¹⁴  AND  |Ξ · u_m| > 10⁻¹⁰:  the state is floppy → return g = 0
```

### What each implementation does

| Function | Threshold | Floppy check |
|---|---|---|
| `calculateShearModulusWithEigenDecomposition` (single-prec, ~line 2898) | **none** — all modes included | **no** |
| `calculateShearModulusWithEigenDecompositiond` (double-prec, ~line 3116) | `> 1e-5` | **no** |
| Python `script_TO_CALC_SM_FROM_CPP.py` | `> 1e-4` | **no** |

### Why this matters

- The single-precision function includes every mode above the LU-rank cutoff with no
  threshold, so near-zero eigenvalues divide `(Ξ · u_m)²` by a tiny number and
  dominate the nonaffine correction.
- Using `1e-5`, `1e-4`, and no threshold produces three different values for the same
  configuration, especially near the rigidity transition where small eigenvalues matter.
- The document's floppy-state check (detect states where g should be exactly zero) is
  absent everywhere.

### Proposed fix

1. **C++ double-precision function:** keep `1e-5` or refine to `1e-10` given double
   precision; add floppy-state early-return:
   ```cpp
   if (currenteigenvalue > 1e-10) {
       nonaffine_shearportion += ...;
   } else if (fabs(xi_dot_u) > 1e-7) {
       return 0.0;  // floppy state
   }
   ```
2. **C++ single-precision function:** add an equivalent `if (currenteigenvalue > 1e-5)`
   guard before the contribution is accumulated (it currently has none).
3. **Python script:** change `1e-4` to match the C++ threshold (`1e-5` or `1e-10`).

---

## Issue 4 — `CalculateShearModulusForPythonCalculation` exports a float (32-bit) Hessian

### What the code does (~line 3180)

```cpp
MatrixXf Hessian(2 * numofvertices, 2 * numofvertices);   // float32 !
Hessian = calculateHessian(pythonfile);
// ... then writes Hessian to pythonfile
```

Python reads this float32 matrix and calls `np.linalg.eigh` (float64 internally), but
every element already has only ~7 significant digits of precision.

### Why this matters

Near the rigidity transition the relevant eigenvalues can be O(10⁻⁶) to O(10⁻⁸). With
a float32 Hessian the relative rounding error in a matrix element is ~10⁻⁷, so elements
with magnitude ~1 already corrupt eigenvalues at the 10⁻⁷ level. The very eigenvalues
that control the non-affine correction are then noise-dominated.

### Proposed fix

Switch to the double-precision Hessian in this function:

```cpp
// In CalculateShearModulusForPythonCalculation:
MatrixXd Hessian(2 * numofvertices, 2 * numofvertices);   // double
Hessian = calculateHessiand(pythonfile);                  // double version
```

---

## Issue 5 — Normalization factor inconsistency between the two C++ shear-modulus functions

### What each function uses

**Single-precision `calculateShearModulusWithEigenDecomposition` (line 2927):**
```cpp
double normfactor = 1 / (K * pow((*this).A0, 2) * (*this).cells.size());
```

**Double-precision `calculateShearModulusWithEigenDecompositiond` (line 3145):**
```cpp
double normfactor = (*this).cells.size() / (K * pow(LX * LY, 2));
```

### Why they differ

Let `Nc = cells.size()` and `A_box = LX*LY`.

- Single: `1 / (K · A0_tissue² · Nc)`
- Double: `Nc / (K · A_box²)`

These are equal only when `(*this).A0 == A_box / Nc`, i.e., the tissue-level member
`A0` equals the average per-cell box area. If `(*this).A0` stores the **target** area
(which may differ from the actual average area, e.g., after cell division or anisotropic
growth), the two functions return different shear moduli for the same configuration.

### Proposed fix

Pick one canonical normalization and use it everywhere.  The standard dimensionless
shear modulus in the vertex-model literature is `g* = g / (K · A₀)`.  Using the
actual box area to define A₀:

```cpp
// Consistent normalization: g* = (Nc / (K * A_box²)) * (d²E/dγ² - nonaffine)
double A_box = LX * LY;
double normfactor = (*this).cells.size() / (K * A_box * A_box);
```

Apply this same expression to both functions and to `CalculateShearModulusForPythonCalculation`.

---

## Issue 6 — `dA/dγ` function naming is inverted (maintenance hazard)

### What the code does (lines 2838/3052)

```cpp
dadgcheck    = (*agit).GetdAdgammaForSimpleShearAlongXAxisALTERNATE();  // USED in calc
dadgcheckALT = (*agit).GetdAdgammaForSimpleShearAlongXAxis();            // logged, unused
```

The function named "ALTERNATE" is the one used for all downstream calculations
(`dadgterm`, `BoundaryVector`). The result of the function without the "ALTERNATE"
suffix is logged but never used.

### Why this matters

If the two functions return different values (possible because they use different
algorithms — one integrates along crossing edges, the other computes vertex averages),
the wrong one may be active. The swapped naming obscures which algorithm is authoritative.

### Proposed fix

1. Determine which function correctly implements `dA_α/dγ = L_y · [ȳ_L^phys − ȳ_R^phys]`
   (document §4) by running a finite-difference check: perturb γ by ε, recompute A,
   and compare (A(γ+ε) − A(γ−ε)) / (2ε) to both functions.
2. Rename the correct function to `GetdAdgammaForSimpleShearAlongXAxis` (removing
   ALTERNATE) and delete or clearly mark the backup as deprecated.

---

## Issue 7 — `d²P/(dγ · dr_i)` mixed derivative: sign convention needs finite-difference verification

### What the document says (§6, Xi formula)

```
Ξ_i^μ = Σ_{α ∈ cells(i)} [  K_A·(dA_α/dγ)·(dA_α/dr^μ_i)
                             + K_P·(dP_α/dγ)·(dP_α/dr^μ_i)
                             − Π_α · d²A_α/(dγ dr^μ_i)
                             + T_α  · d²P_α/(dγ dr^μ_i) ]
```

### What the code does (lines 3086–3088)

```cpp
d2Pdrdgamma = (*agit).Getd2PdgammadvForSimpleShearAlongXAxis((*vit));
BoundaryVector(2*vmapid,   0) += tension_magnitude * d2Pdrdgamma.x;
BoundaryVector(2*vmapid+1, 0) += tension_magnitude * d2Pdrdgamma.y;
```

The area-stress term (`cell_pressure * d2Adrdgamma`) has a verified correct sign because
`cell_pressure = −Π` (line 3045). The perimeter-stress term depends on the sign and
magnitude returned by `Getd2PdgammadvForSimpleShearAlongXAxis`, which contains non-trivial
LEBC sign conventions (CCWVertFactor, updownfactor, LY/l factors).

### Proposed verification

Compute Xi numerically:
```
Ξ_i^μ ≈ [dE(γ+ε)/dr^μ_i − dE(γ−ε)/dr^μ_i] / (2ε)  evaluated at the equilibrium config
```
Compare to the assembled `FullBoundaryVector`. A mismatch would reveal a sign or
magnitude error in `Getd2PdgammadvForSimpleShearAlongXAxis`.

---

## Issue 8 — `d2Pdridrj_ALT` (edge-based perimeter-stress Hessian) is computed but never used

### What the code does (lines 2688/2322 in `calculateHessiand` / `calculateHessian`)

```cpp
d2Pdridrj_ALT = (*this).calculated2PdridrjMatrixd(errorfile);   // computed …
Hessian = d2Adridrj + d2Pdridrj + dAdridAdrj + dPdridPdrj;     // … but not used
```

The edge-based function `calculated2PdridrjMatrixd` computes the perimeter stress
Hessian using `edge_tension_mag` (the per-edge tension, summed from both adjacent cells)
and the rotated-edge-vector formula that is algebraically equivalent to the document's
`−T/l · (I − t̂⊗t̂)`. The cell-by-cell `d2Pdridrj` uses a separate per-cell accumulation.

The commented-out debug line `"d2Pdridrj alt minus standard"` shows these were compared
during development but the comparison result was never documented.

### Why this matters

If the two computations disagree on shared edges (which have contributions from two
cells), one of them is wrong. Using the wrong one silently corrupts the Hessian, which
in turn corrupts the shear modulus.

### Proposed fix

1. **Uncomment** the debug line in `calculateHessiand`:
   ```cpp
   *errorfile << "d2Pdridrj alt minus standard = \n"
              << d2Pdridrj_ALT - d2Pdridrj << "\n" << endl;
   ```
   Run on a small test case and confirm the difference is zero (or machine-epsilon).

2. If they agree, remove the ALT computation to avoid dead code.

3. If they disagree, determine which is correct via a finite-difference check
   (`H_ij ≈ (f_i(r_j+ε) − f_i(r_j−ε)) / (2ε)`) and switch to the correct one.

---

## Summary Table

| # | Location | Severity | Nature |
|---|----------|----------|--------|
| 1 | `vertex_model_master_derivations.tex` §4 | High | Document formula wrong for LEBC (`dP/dγ`) |
| 2 | `vertex_model_master_derivations.tex` §4 | High | Document formula wrong for LEBC (`d²P/dγ²`) |
| 3 | `AGPerimTissue.cpp` lines 2898/3116; Python script | High | Inconsistent eigenvalue thresholds; missing floppy-state check |
| 4 | `AGPerimTissue.cpp` line 3180 (`CalculateShearModulusForPythonCalculation`) | High | Exports float32 Hessian — precision loss near transition |
| 5 | `AGPerimTissue.cpp` lines 2927 vs 3145 | Medium | Normalization differs between single- and double-precision functions |
| 6 | `AGPerimTissue.cpp` lines 2838/3052 | Medium | `dA/dγ` function names inverted; wrong function may be active |
| 7 | `AGPerimTissue.cpp` lines 3086–3088 | Medium | `d²P/(dγ dr_i)` signs unverified by finite difference |
| 8 | `AGPerimTissue.cpp` lines 2688/2322 | Low | Edge-based Hessian ALT computed but silently discarded |
