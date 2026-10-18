# Vertex-model code changes for the Hessian shear modulus: a detailed report

**Branch:** `MondalShearDev` · **Period:** 2026-10-10 to 2026-10-15 · **Status:** library changes committed (not pushed); the quartic re-equilibration is still running (see §8).

This report records every change made to the vertex-model C++ code for the Hessian shear-modulus calculation. For each change it gives the symptom, the root cause, the code before and after, how the change was verified, and its effect on results. It also defines how floppy (fluid) tissues are handled and reports the measured rate of floppiness. It is intended as the methods and erratum record for publication.

---

## 0. Summary

| # | What | Where | Commit | Effect on G |
|---|---|---|---|---|
| 1 | Unified the eigenvalue threshold and added a floppy-state return (float + double); fixed the float-truncated eigenvalue; unified the normalization | `AGPerimTissue.cpp`, `calculateShearModulusWithEigenDecomposition[d]` | `5f7a19504` | Makes the float and double paths consistent; normalization uses the true box area |
| 2 | **Sign of the area-stress block of the Hessian** (main bug) | `AGPerimTissue.cpp`, `calculateHessian[d]` | `f09065e4b` | Removes spurious soft modes, false floppy returns and negative G near the transition |
| 3 | Lees–Edwards dA/dγ: wrong sign for ~2.7% of boundary-crossing cells | `AGCell.cpp`, `GetdAdgammaForSimpleShearAlongXAxisALTERNATE` | `7a398ed78` | Corrects Ξ and the affine term |
| 4 | Lees–Edwards d²A/dγ dy: wrong sign for every corner cell | `AGCell.cpp`, `Getd2AdgammadyForSimpleShearAlongXAxis` | `7a398ed78` | Corrects Ξ |
| 5 | **C4 shadowing**: forces always used C4 = 1 | `AGTissue.h` | `4090de861` | Quartic tissues with C4 ≠ 1 must be re-equilibrated (§5, §8) |
| 6 | Zero-mode and floppy thresholds made relative to the matrix scale; every mode examined | `AGPerimTissue.cpp`, `calculateShearModulusWithEigenDecomposition[d]` | `9b39943bb` | Removes random, round-off-driven floppy flags. Solid-side G is unchanged |
| D1 | New batch driver for the Hessian G | `AGPerimHessian_test/AGPerimHessian_test/MainHess_claude.cpp` | `143e54893` | No physics change |
| D2 | New batch driver for equilibration | `AGPerimVertex/MainAGPerim_claude.cpp` | `e52d1279c` (+ `_NE.txt` output 2026-10-16, chunked restarts and save-on-signal 2026-10-17) | No physics change |

**Lesson for future code:** every bug above came from writing the code like a derivation, i.e. assuming exact arithmetic and conventions that were never tested numerically. See §11.

**End result.**
- After changes 1–4, the library's G matches an independent finite-difference (FD) reference to ~5e-9 relative for every T2ss tissue with p0 ≤ 3.82 (128 tissues).
- Across 22 data sets (low/high polydispersity × quadratic/log/exp potentials × A0 = 0.5, 1, 2 ⟨A⟩, plus SWAB and the C4 = 1 quartic sets), the Hessian G matches the squeeze-and-fit benchmark: G_Hessian / G_squeeze = 1.00–1.04 in the solid (p0 ≤ 3.80).

**Line-ending note.** `git show --stat` overstates commits 1 and 5 (7,073 and 298 lines). That is CRLF/LF churn: the working tree has pre-existing line-ending drift. The real edits, from `git show -w --ignore-cr-at-eol`, are:

| Commit | Real lines added | Real lines removed |
|---|---|---|
| 1 | 50 | 7 |
| 2 | 20 | 8 |
| 3 + 4 | 98 | 198 |
| 5 | 7 | 1 |
| 6 | 34 | 6 |

---

## 1. Background: what the code computes

The shear modulus is computed by linear response:

  G = (1/A_box) [ ∂²E/∂γ² − Ξᵀ H⁺ Ξ ],  Ξ_i = ∂²E/∂r_i∂γ,  H_ij = ∂²E/∂r_i∂r_j

- γ is a Lees–Edwards (LEBC) simple shear along x. The periodic image above the box shifts by γL_y.
- H⁺ is the pseudo-inverse: modes with near-zero eigenvalue are dropped.
- The library returns `Smd = G · N_c / (K (L_X L_Y)²)`. Every result here is converted back with G = Smd · K L_X L_Y / N_c.

Per cell, the energy is E_c = ½ K A0 f(A/A0) + ½ Λ (P − P0)². The area potential f is chosen per cell by `pressureType`:

| `pressureType` | f(a), a = A/A0 |
|---|---|
| quadratic (default) | (a−1)² |
| `LOG` | −2 log(1 + a) |
| `EXPONENTIAL` | 2 exp(−a) |
| `QUARTIC` | (a−1)² + C4 (a−1)⁴ |
| `QUARTIC2.0` | (1−6C4)(a−1)² + C4(a−1)⁴ |

The pressure Π = −∂E/∂A is computed by `AGCell::GetPressure`, and the area stiffness ∂²E/∂A² by `AGCell::GetNegativeDPressureDArea`.

**Verification tool used throughout.** An independent re-implementation of the energy, frozen at each tissue's topology. All derivatives (H, Ξ, ∂²E/∂γ²) are computed from it by 4-point central finite differences. The code is in `claude_VM_runs/T2ss_SM_Hess_testing/src/{MainHessFDCheck.cpp, RefModel.h}`. It was extended with the quartic branch in `Quartic_C4fix/src/`.

---

## 2. Change 1: eigenvalue threshold, floppy state, float truncation, normalization (`5f7a19504`, 2026-10-10)

`AGPerimTissue.cpp`, `calculateShearModulusWithEigenDecomposition` (float) and `...Decompositiond` (double).

**Problems**
- (a) The float version had **no** eigenvalue threshold: every mode inside the LU-rank range entered as ξ²/λ, including near-zero ones.
- (b) The double version stored the eigenvalue in a `float` before comparing it to the threshold.
- (c) The two versions used different normalizations: `1/(K A0² N_c)` versus `N_c/(K (L_X L_Y)²)`. They agree only if A0 = L_X L_Y / N_c.

**Before** (float version)
```cpp
nonaffine_shearcontribution = 1 / hessianeigen.eigenvalues()(index, 0) * pow((u_index^T * Xi), 2);   // every mode, no cut
...
double normfactor =  1 / (K * pow((*this).A0, 2) * (*this).cells.size());
```
**Before** (double version)
```cpp
float currenteigenvalue;            // truncated before the comparison
if (currenteigenvalue > 1e-5) { ... }
```
**After** (both versions)
```cpp
const double EIGENVALUE_THRESHOLD = 1e-5;   // later made relative, see Change 6
const double FLOPPY_TOLERANCE = 1e-4;
double xi_dot_u = (u_index^T * FullBoundaryVector)(0,0);
if (currenteigenvalue > EIGENVALUE_THRESHOLD)      nonaffine += xi_dot_u*xi_dot_u/currenteigenvalue;
else if (fabs(xi_dot_u) > FLOPPY_TOLERANCE)        return 0.0;      // "floppy state": G = 0
...
double normfactor = (*this).cells.size() / (K * pow(LX * LY, 2));   // both versions
```
and `double currenteigenvalue;` in the double version.

**Verification.** Full build and link. The FD check later confirmed the double-precision path (§3).

---

## 3. Change 2: sign of the area-stress block of the Hessian (`f09065e4b`, 2026-10-14). Main bug

`AGPerimTissue.cpp`, `calculateHessian` (float) and `calculateHessiand` (double).

**Symptom.** Spurious soft modes at p0 ≈ 3.80–3.86. They produced false floppy returns and negative G near the rigidity transition.

**Root cause**
- The Hessian contains the pressure ("stress") term −Π_α ∇_i∇_j A_α.
- For the positive polygon area A = ½ Σ_k (x_k y_{k+1} − x_{k+1} y_k) with k → k+1 **counterclockwise**, ∇_i∇_{i+1} A = ½[[0, 1], [−1, 0]].
- The code loops over `SavedInfoForEvolve_CWVerts`, which is **clockwise**, but used the counterclockwise signs (the code comment said "SIGNS SET BY CCW CONVENTION"). Every off-diagonal area-stress block therefore had the wrong sign.
- Because the (i, j) and (j, i) blocks are transposes, the Hessian stayed symmetric, so a symmetry check could not catch this.
- It vanishes when Π = 0 and matters most near the transition, where the perimeter stiffness is small.

**Before**
```cpp
d2Adridrj_percell(2*jm1_id,   2*j_id+1) += 0.5*cell_pressure;
d2Adridrj_percell(2*jm1_id+1, 2*j_id  ) -= 0.5*cell_pressure;
d2Adridrj_percell(2*jp1_id,   2*j_id+1) -= 0.5*cell_pressure;
d2Adridrj_percell(2*jp1_id+1, 2*j_id  ) += 0.5*cell_pressure;   // "SIGNS SET BY CCW CONVENTION"
```
**After** (all four signs flipped, for the clockwise list)
```cpp
d2Adridrj_percell(2*jm1_id,   2*j_id+1) -= 0.5*cell_pressure;
d2Adridrj_percell(2*jm1_id+1, 2*j_id  ) += 0.5*cell_pressure;
d2Adridrj_percell(2*jp1_id,   2*j_id+1) += 0.5*cell_pressure;
d2Adridrj_percell(2*jp1_id+1, 2*j_id  ) -= 0.5*cell_pressure;   // signs for the CW vertex list
```

**Verification**
- Before the fix, the relative error of this block against the FD Hessian was exactly 2, the signature of a pure sign flip.
- After it, the full Hessian matches FD to 1.8e-11 (median over 128 T2ss tissues).
- The spurious soft modes and false floppy returns at p0 3.80–3.82 are gone.

**Origin.** The derivation notes (`vertex_model_master_derivations`, area-stress term) give the counterclockwise formula without saying so. A correction is written into the notes.

---

## 4. Changes 3 and 4: Lees–Edwards area derivatives (`7a398ed78`, 2026-10-14)

`AGCell.cpp`. A new static helper `HelperLEBCUnwrap` walks the clockwise vertex list with `Edge::GetVectorDifference`. It returns three things:
- the unwrapped polygon;
- each vertex's periodic image level ℓ_k, the running sum of the box heights n_y crossed by each edge;
- σ = ±1, which makes the shoelace area positive.

Under LEBC, a vertex at image level ℓ_k moves by ℓ_k γ L_y in x.

### Change 3: dA/dγ (`GetdAdgammaForSimpleShearAlongXAxisALTERNATE`)

**Before**
- Took the y-midpoints y₁ and y₂ of the first two edges that cross the top/bottom boundary.
- Returned L_y(y₁ − y₂) or L_y(y₂ − y₁), choosing the sign from `crosscellvec.x`, the x component of a partial edge sum.
- The magnitude was right. The sign heuristic failed for **41 of 1546** boundary-crossing cells (2.7%), including cells that are not corner cells.
- It also assumed exactly two crossing edges.

**After**
```cpp
dA/dγ = L_y Σ_k ℓ_k ∂A/∂x_k ,   ∂A/∂x_k = σ · ½ (ỹ_{k+1} − ỹ_{k−1})   (unwrapped coordinates)
```
This is the exact chain rule of the shoelace area. It is 0 for cells that do not cross the boundary, and it handles corner cells and multiple crossings. For a cell that crosses once, it reduces to L_y(ȳ_L − ȳ_R), where R and L are the crossing edges walked upward and downward in counterclockwise order.

### Change 4: d²A/dγ dy_j (`Getd2AdgammadyForSimpleShearAlongXAxis`)

**Before.** Wrong sign for **every** corner cell (cells crossing both the x and y boundaries), 122 of 122.

**After**
```cpp
∂²A/∂γ∂y_j = σ (L_y/2) (ℓ_{j−1} − ℓ_{j+1}),   ∂²A/∂γ∂x_j = 0
```

**Verification**
- 0 bad cells in 128 tissues.
- Ξ matches FD to 1e-10.
- After Changes 1–4, G_code = G_FD to ~5e-9 relative for every tissue with p0 ≤ 3.82.

**Note on conventions.** Under LEBC, ∂A_α/∂γ ≠ 0 for each boundary-crossing cell; only the sum over cells vanishes. So Ξ contains a K_A ∇A ∂A/∂γ term and ∂²E/∂γ² contains a K_A(∂A/∂γ)² term, and the code includes both. A homogeneous affine shear r → Λr instead gives ∂A_α/∂γ = 0 for every cell. The two conventions give the same G, but their intermediate quantities differ. The `.md` derivation notes mix the two conventions, and corrections are written there.

---

## 5. Change 5: the C4 shadowing bug (`4090de861`, 2026-10-15)

`AGVertexModel/AGTissue.h`.

**Symptom, found before any quartic run**
- (i) In every quartic EQ file, the stored per-cell pressures imply C4 = 1.000, whatever C4 the file header gives.
- (ii) Data sets that differ only in C4 have **identical vertex positions**:
  - QuTw (C4 = 3), Qu2w (2), QuOw (1) and MTQ3 (10);
  - QuTs (30), Aq2w (2) and AqOw (1);
  - Lq2w (2) and LqOw (1).
- (iii) The stored `finalEnergy` matches the header C4.

So the forces never saw C4, while the energy did.

**Root cause.** C4 was declared twice:
```cpp
// AGConstants.h
extern double C4;                 // global, defined in every Main*.cpp as  double C4 = 1;
// AGTissue.h  (class AGTissue)
double C4;                        // MEMBER with the same name
```
Inside tissue member functions, the bare name `C4` resolves to the member:
- `AGPerimTissue::ReadInputData` set the **member** from the file;
- `CalculateFreeEnergy` read the **member**;
- `AGCell::GetPressure` and `GetNegativeDPressureDArea` are not tissue members, so they read the **global**. Nothing ever assigned it, so it kept the `Main*.cpp` default of 1. These two functions drive the forces during equilibration and squeezing, and the Hessian.

**Before → after**
```cpp
-	double C4; // for the quartic potential
+	// double C4; ... <- REMOVED 2026-10-15 (C4 shadowing bug). [explanatory comment]
```
With the member removed, `ReadInputData` sets the single global C4, and every function uses it. No code referred to `tissue.C4` explicitly.

**Why this pattern breaks and others don't.** Parameters that work are either:
- true per-cell members, filled from the `Cells` block and used by cell methods (A0, P0, `pressureType`); or
- true globals with a single declaration (K, GAMMA, LAMBDA).

C4 was both, so the reader and the user of the value saw different variables. `equilibriumtolerance` and `energyListSize` also exist as both a `Main*.cpp` global and an `AGPerimTissue` member. They were checked: they are read and used through the same member, so they are consistent.

**Verification** (`Quartic_C4fix/`)
- **FD check with an independent quartic reference energy,** E_A = ½(K/A0)(A−A0)² + ½(K/A0)C4(A−A0)⁴/A0², written from the energy alone and not from the library's pressure functions:

| | C4 = 1 (QuOw) | C4 = 3 (QuTw) | C4 = 30 (QuTs) |
|---|---|---|---|
| library `GetPressure` vs −dE/dA | 3e-16 | 4e-16 | 9e-16 |
| Hessian vs FD | 1.8e-11 | 1.8e-11 | 1.8e-11 |
| Ξ vs FD | 5–9e-11 | 1.0e-10 | 6e-11 |
| G_code vs G_FD | 2e-9 | 3e-9 | 1e-9 |
| net force on the OLD EQ file under the true C4 | 2e-11 | **10.2** | **1.75** |

  The last row quantifies the old bug: tissues with C4 ≠ 1 were far from force balance under their own energy.
- **Regression.** Re-equilibrating QuOw tissues (C4 = 1, where the bug changes nothing) from the raw inputs with the current repo code reproduces the 2023 Great Lakes EQ files: final energy to 1e-13, vertex positions to 1e-5. This shows the repo library evolves exactly like the 2023 code.
- **Non-quartic results are unchanged.** The T2ss check tissue gives an identical `Smd`.
- **Pilots.** In newly equilibrated tissues, the pressures imply exactly the header C4 (3, 10, 30).

**Consequences for existing data**
- The C4 ≠ 1 quartic sets (QuTs, QuTw, Qu2w, Aq2w, Lq2w, MTQ3) were equilibrated with C4 = 1 forces. Their squeeze benchmarks relaxed with C4 = 1 forces but recorded the energy with the file's C4, so **those squeeze moduli are not the modulus of either model** (e.g. QuTw's G is non-monotonic in p0).
- The C4 = 1 sets (QuOw, AqOw, LqOw) are valid.
- The C4 ≠ 1 sets are being re-equilibrated (§8). By the user's decision, quartic G is computed by the Hessian method only; no new squeeze runs.

---

## 6. Change 6: relative zero-mode and floppy thresholds (`9b39943bb`, 2026-10-15)

`AGPerimTissue.cpp`, both shear-modulus functions.

**Symptom.** Physically identical tissues were classified differently.
- In the quadratic model with K ∝ A0, at fixed box area, A0 drops out of the energy because Σ A_α is fixed. So L2Ts, M2Ts and T2ss share one Hessian and one Ξ (both checked: λ_max and ‖Ξ‖ agree).
- Yet the library flagged 38 T2ss tissues as floppy against 47 for M2Ts.
- Some deep-fluid tissues returned "non-floppy" with G/K ~ 1e-8.

**Root cause** (two parts)
1. **The loop skipped part of the spectrum.** It started at `index = size − lu_decomp.rank()`, so every mode that Eigen's `FullPivLU` counted as null space was **skipped without a floppy test**. A fluid tissue has ~N_v − 2N_c ≈ 200 genuine zero modes. Which of them fell inside the examined range depended on round-off in the LU rank.
2. **Absolute cuts at the noise level.** The cuts were 1e-5 for λ and 1e-4 for ξ·u. Large uniform pressures (A far from A0 in L2/M2) leave round-off of ~1e-8·λ_max in the near-zero eigenvalues. The M2Ts "floppy" trigger was ξ·u = 4e-4 ≈ **4e-9 ‖Ξ‖**, i.e. noise.

**Physics.** In a fluid (zero-energy) state, Ξ lies in the range of H, so the nonaffine term cancels the affine term exactly and G = 0. A genuine floppy coupling of the shear to a zero mode is a different situation, and the noise-triggered flag mixed the two.

**Before**
```cpp
const double EIGENVALUE_THRESHOLD = 1e-5;
const double FLOPPY_TOLERANCE = 1e-4;
for (int index = size - lu_decomp.rank(); index < size; index++)
```
**After**
```cpp
const double REL_EIGENVALUE_THRESHOLD = 1e-8;
const double REL_FLOPPY_TOLERANCE    = 1e-6;
const double EIGENVALUE_THRESHOLD = REL_EIGENVALUE_THRESHOLD * max|lambda|;
const double FLOPPY_TOLERANCE    = REL_FLOPPY_TOLERANCE    * |Xi|;
for (int index = 0; index < size; index++)          // every mode is examined
```
- λ ≤ cut: zero mode, dropped from the pseudo-inverse.
- λ ≤ cut **and** |ξ·u| > 10⁻⁶‖Ξ‖: floppy return, G = 0.
- Otherwise the mode enters as ξ²/λ.
- With λ_max ≈ 630 for these tissues, the relative eigenvalue cut equals the old 1e-5, so verified results do not move. The floppy cut is ~250× above the measured round-off.

**Verification** (2,446 tissues re-run, `Quartic_C4fix/relthr/`)
- G for every tissue with p0 ≤ 3.80: **unchanged** (maximum relative change 0).
- T2ss at p0 3.81–3.83: G unchanged and equal to G_FD. 20 of 21 tissues agree to all printed digits; one at 3.83 sits at the 3e-4 noise level.
- Library floppy returns: 776 → 25.
- L2Ts / M2Ts / T2ss now classify every tissue identically (0 disagreements out of 112).
- Geometrically identical tissue pairs agree in G to 4e-9 (median) and 2.6e-5 (max).

---

## 7. Floppiness: definition, averaging, measured rate

**Averaging rule** (user decision, 2026-10-15)
- A tissue whose library call returns the floppy state is a fluid with **G = 0**, and it enters the seed average as 0.
- Every other tissue enters with its computed G.
- No tissue is dropped unless it has no result at all. Every (data set, p0) point therefore has a value; the summary tables were checked and contain no blanks.
- Mean ± SEM over the 7 seeds (rng 100–160), SEM = std(ddof = 1)/√n.
- Plots show G/K_plot with K_plot = 50/√1250, as in the MATLAB scripts.
- The squeeze benchmark keeps its own MATLAB convention, in which 0 means missing.
- Superseded: the MATLAB convention of dropping zeros had removed whole fluid-side points (e.g. LqOw at p0 3.89/3.90) for the Hessian. The non-floppy-only mean is still tabulated separately (`G_hess_nonfloppy_*`).

**Rate of floppiness**
- After Change 6, fluid tissues mostly come back from the pseudo-inverse with small but **reproducible** residuals, G/K ~ 1e-8 to 1e-3, rather than an exact zero. The residual comes from soft modes just below the eigenvalue cut. It is not noise: identical tissues reproduce it to ~1e-8.
- There is therefore no natural gap, and calling such a tissue fluid requires a cut.
- Main definition: **fluid if floppy or |G|/K < 1e-4**. The fractions for cuts of 1e-5 and 1e-3 are reported as a sensitivity band.
- Changing the cut by a factor of 10 moves the 50%-fluid point by ~0.01 in p0.

| Set (first p0 with ≥ 50% fluid, cut 1e-4) | |
|---|---|
| Low polydispersity, all potentials and A0 (T2ss, L2Ts, M2Ts, LoTs, L2Lo, M2Lo, ExTs, L2Ex, M2Ex, SWAB) | 3.84–3.85 |
| Quartic C4 = 1 (QuOw, AqOw, LqOw) | 3.85 |
| High polydispersity (`*_LP`) | 3.88 |

- Fluid tissues first appear at 3.82–3.83 (low polydispersity) and 3.84–3.86 (high).
- The full per-p0 table, with both cuts and the non-zero-seed mean G, is `Hess_SM_groups/analysis/floppy_fraction.tsv`. The figure is `floppy_fraction.png`.

**Caveat.** On the fluid side, only a few seeds per p0 have non-zero G, so the rate and those G values rest on few tissues. More seeds are planned (todo list).

---

## 8. Drivers and data handling

**D1. `MainHess_claude.cpp`** (copy of `MainHessianTest.cpp`, commit `143e54893`)
- argv input and an output folder.
- **No evolution**: the Hessian is taken at the EQ state the squeeze benchmark started from.
- Sets each cell's `SavedInfoForEvolve_CWVerts`. `Evolve` normally does this, and without it G was ~96× wrong in the smoke test.
- `Smd` calls the double-precision function (the original called the float one twice).
- Debug dumps and the Python-only call removed.
- Writes `<tissue>_SM.txt` (line 1 Smd, line 2 Smf).

**D2. `MainAGPerim_claude.cpp`** (copy of `MainAGPerim.cpp`)
- argv input and an output folder.
- No per-step `_out.txt` dump.
- Writes `_EQ.txt` only if the tissue equilibrated, with the initial/final energy footer, as the 2023 Great Lakes driver did.
- (2026-10-16) If the tissue did not equilibrate by TOTALTIME, it writes the final state to `_NE.txt`, as the user's Great Lakes `MainAGDebug.cpp` does. That file is a complete input file, so a run can be continued from it. Checked: the restart's initial energy equals the saved final energy exactly, and a rerun reproduces the original run's log bit for bit.
- (2026-10-17) Chunked, restartable runs (protocol: `vm_claude_plans/2026-10-17_chunked-equilibration-protocol.md`):
  - An optional third argument sets the chunk length, overriding the file's TOTALTIME. The batch file `batch_EQ_chunked.sbat` runs chunks of 1000 from the latest `_NE.txt` until the tissue equilibrates or reaches a total of t = 10000.
  - CURRENTTIME in the output is the tissue's total evolved time (start time read from the input + time run in this chunk).
  - SIGUSR1/SIGTERM set a flag that ends the loop after the current `Evolve` call, and the state is written to `_NE.txt`. SLURM sends USR1 15 min before the time limit, so a run that hits the limit still saves its state. Tested in a 10-min SLURM job.
  - Output names stay `<tissue>_ELS_3_EQ.txt` / `_NE.txt` over any number of chunks (the previous `_NE.txt` is kept as `_NE_prev.txt`). `_error.txt` is appended with a header per chunk. Each chunk adds one line to `<tissue>_ELS_3_chunks.tsv` (start/end time, energies, final slope, outcome, wall time).
  - Each chunk restarts the energy-slope window and the step size, as any restart from a file does. Checked: the restart's initial energy equals the saved final energy exactly.
- (2026-10-17) 36 tissues still running under the old binary were stopped with gdb (`set var tissue.equilibrated = 1`, user-approved), so the old driver wrote its state. These files were renamed to `_NE.txt` with the footer corrected and continued under the chunked protocol. The originals are kept as `_gdbstop_EQ.txt`.
- The evolution loop and the equilibrium criteria are unchanged.
- Observation, not changed: in this re-equilibration, 544 of 607 tissues met the slope test (< 1e-15) only after the step size collapsed and the last 3 energies became identical, so the regression slope was exactly 0. The slopes just before were 1e-13 to 1e-9. This is the existing termination behavior of `AGPerimTissue::Evolve`.

**Data overwrite** (user-approved)
- The old C4-bug EQ files of QuTs, QuTw, Qu2w, Aq2w, Lq2w and MTQ3 (652 files) are kept in `CompressionInducedFluidizationData/claude_hessian_testing/Quartic_C4fix/original_EQ_C4bug/`.
- After re-equilibration (SLURM job 64105735), the new EQ files replace them in place.
- Old EQ files of tissues that fail to equilibrate are moved to that folder.
- Squeeze files and `_dE` folders are left untouched.

---

## 9. Checked and found correct (no change)
- **Perimeter terms of the Hessian.** Including d²P/dr dr: the standard and `_ALT` versions agree to 1e-12, so the ALT copy is redundant.
- **Perimeter strain derivatives.** ∂P/∂γ, ∂²P/∂γ² and ∂²P/∂γ∂r are correct; the sign convention is s = +1.
- **Pressure and stiffness functions.** For all potentials, including both quartic variants, the pressure is the exact derivative of the energy and the stiffness the exact derivative of the pressure.
- **The ~2% offset between Hessian G and squeeze-fit G in the low-polydispersity sets is not a code error.** It comes from each finite sample's pure-vs-simple shear anisotropy plus a ~1% bias from fitting only on one side (+ε). The benchmark's energies are reproduced to 6 digits.

## 10. Where things are
- **Commits:** `5f7a19504`, `f09065e4b`, `7a398ed78`, `4090de861`, `9b39943bb` (library); `143e54893` (D1).
- **Work folders:** `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/{T2ss_SM_Hess_testing, Hess_SM_groups, Quartic_C4fix}/` (each has a README), backed up to `/nfs/turbo/lsa-dkluben/abmondal/CompressionInducedFluidizationData/claude_hessian_testing/`.
- **Derivation notes:** `vertex_model_master_derivations.{md,tex}` carry ALL-CAPS "CLAUDE CORRECTION" notes for the author to rewrite.

## 11. Lesson for future code: don't code like a theorist

Every bug in this report has the same root cause. The code was written as a transcription of a derivation, with exact arithmetic and the derivation's conventions taken for granted. It was never checked numerically against the energy it was supposed to differentiate. Each bug was locally plausible and invisible on paper. Once a finite-difference check existed, each one showed up directly as a mismatch.

| Theorist's assumption | Where it broke | What to do instead |
|---|---|---|
| "The derivation uses counterclockwise order, so the code does too" | Change 2: the d²A/dr dr signs were written for a counterclockwise list, but `cwverts` is clockwise. The comment even said "SIGNS SET BY CCW CONVENTION". | State each convention as a property of the data structure and check it in code (e.g. the sign of the shoelace sum). Never infer it from the notes. |
| "Zero is zero" | Changes 1 and 6: absolute cuts (1e-5, 1e-4) sat at round-off in some tissues and above it in others. `FullPivLU`'s rank, used to skip null modes, depended on noise, so identical Hessians were classified differently. | Make every tolerance relative to the scale of the quantity (λ_max, ‖Ξ‖). Never let a decision depend on a numerically computed rank. Log what was cut. |
| "In exact arithmetic, skipping the null space is the pseudo-inverse" | Change 6: the `size − rank` loop start (your 2024 code) is correct in exact arithmetic. In floating point it silently skipped modes, including any negative ones and any zero mode that couples to the shear, so they never reached the floppy check. | Examine every mode and decide explicitly, with a stated criterion, what to drop. A mode should never be skipped by index arithmetic. |
| "I can handle the geometry case by case" | Changes 3 and 4: the boundary-crossing area derivatives picked their sign with heuristics (edge-sum direction, the x ordering of intercepts). They failed for 2.7% of crossing cells and for every corner cell. A comment had already flagged "POTENTIAL WORRIES ABOUT CORNER CELLS". | Derive from one general formula (here, the shoelace formula on the unwrapped polygon) instead of special-casing the geometry. A known worry in a comment means a test is due, not a comment. |
| "Single vs double precision doesn't matter" | Change 1: a `float currenteigenvalue` inside the double-precision function cut its eigenvalues down to float precision. | Keep one precision through each pipeline and grep for stray `float` in double code paths. |
| "A parameter is just a parameter" | Change 5: a member `AGTissue::C4` shadowed the global C4. The energy used the file's value and the forces used the default of 1. The two disagreed silently, with no error or warning. | Keep a single source of truth for each physical constant. Compile with `-Wshadow`. Check that forces are the exact negative gradient of the energy for every potential. |
| "The run will finish" | Driver D2 (Claude's own code, 2026-10-15): tissues that did not equilibrate exited without saving their state, so 37 tissues had to be rerun from t = 0. | Every exit path of a long run writes its state. Here that means `_NE.txt`, as the user's `MainAGDebug.cpp` already did. |

**The single most valuable practice: a finite-difference check against the energy.** One independent reference (`MainHessFDCheck.cpp`: energy → FD gradient → FD Hessian → FD Ξ → FD G) exposed changes 2–5 immediately. It also confirmed that the method itself was right, and it is cheap to rerun after any change. For any new energy term or derivative, write the FD check first, then the analytic code. Test it on tissues that hit the edge cases: corner cells, fluid tissues, C4 ≠ 1, every potential.

