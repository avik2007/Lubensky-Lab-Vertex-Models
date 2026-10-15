# Plan: Hessian G for the three low-polydispersity triplets, run from the user's own driver

## Context
The FD check showed that the fixed library's Hessian G matches the squeeze benchmark (T2ss, ExTs). The user now wants the production version of this:
- Run Hessian G with **their own C++ driver**, without the FD/debug machinery.
- Cover the dataset groups their MATLAB overlay script plots together.
- Save the shear moduli.
- Then compare against their saved squeeze results.

The groups come from `/nfs/turbo/lsa-dkluben/abmondal/Avik/PhaseTransitionProject/SM_shearmodulus_fitting_100cells_QuadraticFullPlotScript.m`, lines 1–63. There are three triplets:
- **Quadratic:** L2Ts / M2Ts / T2ss
- **Logarithmic:** L2Lo / M2Lo / LoTs
- **Exponential:** L2Ex / M2Ex / ExTs

The user's decisions:
- Low-polydispersity triplets only.
- **No evolve**: the Hessian is computed on the EQ tissue exactly as loaded.

## Inputs (read-only, nothing modified)
The root is `P = /nfs/turbo/lsa-dkluben/abmondal/Avik/PhaseTransitionProject`, the folder the user's plot script runs from.
- **Unsqueezed EQ files** (`<code>_rng_<seed>_p0_<p0>_ELS_3_EQ.txt`):
  - `P/<code>_EQ/` for T2ss, L2Ts, M2Ts, LoTs, L2Lo and M2Lo.
  - `P/<code>_files/` for ExTs, L2Ex and M2Ex.
  - Each dataset has 7 seeds (100–160) × 16 p0 (3.75–3.90) = 112 files. M2Ex has 111.
- **Squeeze benchmark:** `P/<code>_dE/smdata<code>.mat`. T2ss uses `smdata7.mat`, as the plot script does. The matrices are `p0sm*`, `SMvar*` and `p0byrng*`, covering the same 7 seeds.
- The x axis comes out right: every file has ⟨A⟩ = LX·LY/N = 1250, so the file-name p0 = P0/√⟨A⟩. A0 is 625 for L2, 1250 for T/Lo/Ex and 2500 for M2. K is read from each file's header (0.707 / 1.414 / 2.83).

## Step 1: `MainHess_claude.cpp` (new file in the repo, minimal edits)
- Copy `AG Vertex Model/AGVertexModel/AGPerimHessian_test/AGPerimHessian_test/MainHessianTest.cpp` to `MainHess_claude.cpp` in the same folder. The original is left untouched.
- Only these changes, each with a comment:
  1. **Command-line input:** use `main(int argc, char* argv[])` with `argv[1]` = input file and `argv[2]` = output folder. This replaces the `cin` prompt and the Windows `PathName`. It re-enables the user's own commented-out argv block, which is also the pattern in their Great Lakes copy `~/MainHessianTest.cpp`.
  2. **No evolve:** delete the re-equilibration `while` loop and the extra `Evolve` call (user decision).
  3. **Remove debugging output:** the `_out.txt` / `_EQ.txt` / `_H.txt` / `_P.txt` dumps, the Hessian symmetry printout and the `CalculateShearModulusForPythonCalculation` call. That call serves only the Python path, which is out of scope (lesson L002).
  4. **Fix a bug:** `Smd` currently calls the *float* `calculateShearModulusWithEigenDecomposition` a second time. Change it to `...Decompositiond`.
  5. **New output file:** write `<outdir>/<basename>_SM.txt` with line 1 = Smd (double) and line 2 = Smf (float), at setprecision(20). Line 1 being the modulus matches the user's existing `_SM.txt` convention and their MATLAB reader `ShearModulusfromSMPlotter.m` (`calculateShearModulusfromSMfiles`). The library's verbose diagnostics still go to `_error.txt`, as the original does.
- Library sources are the repo's current, already-fixed `AGVertexModel/*.cpp` and `AGPerimVertex/*.cpp`, used in place and not copied or edited.
- Add the new file to `vm_file_structure.md`.

## Step 2: Build and run (scratch, SLURM)
- **Work folder:** `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Hess_SM_groups/`.
- **Makefile:** model it on the user's `GreatLakes_makefiles/makefile_H`. Same g++ line, with `MainHess_claude.cpp` in place of `MainHessianTest.cpp`, plus `-I` paths to the repo library folders and the vendored Eigen (`AGPerimHessian_test/packages/Eigen3.3.3.9/lib/native/include`). It outputs `agperim_hess_claude`.
  - If `-std=c++11` fails to compile, use `-std=c++14`, the standard my earlier build used, and say so.
- **Parameter file:** `SMparameters_<group>.txt`, one line per tissue: `<EQ path> results/<code>`.
- **Batch files:** `batch_SM_<group>.sbat`, three of them (Ts, Lo, Ex), each an array of about 336 tasks.
  - Modeled on `batch_SM_T2ss.sbat`: `--account=abmondal1`, standard partition, `sed -n ${SLURM_ARRAY_TASK_ID}p`.
  - Settings: `--time=00:30:00`, `--mem=5g`.
- **Smoke test first:** a 2-task array run (1 T2ss and 1 L2Ex tissue). Check that its Smd reproduces the FD-run `G_code_d` for the same T2ss tissue, after the `K·LX·LY/Nc` conversion. Only then submit the full arrays.

## Step 3: Compare and plot
- Script: `analysis/plot_groups_hess_vs_squeeze.py`. It reuses `read_geometry`, `mean_sem` and the MATLAB styling (`K_PLOT = 50/√1250`, fonts) from the earlier `T2ss_SM_Hess_testing/analysis/{compare_to_dE,plot_G_vs_shape_index}.py`, copied over.
- **Conversion:** G_raw = Smd · K·LX·LY/Nc, with K, LX, LY and N from each EQ header. This is the same conversion as the verified driver. A 0 (floppy) result counts as missing, and the number of floppy results per point is reported.
- **Benchmark:** taken straight from `p0sm*` / `SMvar*`, so it is exactly what the user's plots show.
- **Hessian averaging:** mean ± SEM over the same 7 seeds.
- **One figure per triplet:**
  - `G_vs_shapeindex_{Ts,Lo,Ex}_group.{png,pdf}`.
  - The user's colors and markers: L2 black ○, M2 blue □, A0 = ⟨A⟩ red ◇, with the user's legend labels.
  - Squeeze shown as open markers and Hessian as filled markers.
  - The y axis is G/K_PLOT, matching the user's script.
- **Saved values:**
  - `results/<code>/*_SM.txt` (raw, per tissue).
  - `analysis/G_hess_all.tsv`: per tissue, with code, seed, p0, Smd, G_raw and a floppy flag.
  - `analysis/G_groups_summary.tsv`: per dataset and p0, with Hessian mean/SEM/n, squeeze mean/SEM and the ratio.
- **Report:** per dataset, the Hessian/squeeze ratio in the solid regime and where the two diverge.

## Close-out
- Update `vm_claude_actions/vm_claude_todo.md` and `vm_claude_recentactions.md` with a README in the work folder.
- rsync the work folder to `/nfs/turbo/.../claude_hessian_testing/Hess_SM_groups/`, using the same rsync rule as before.
- Commit only `MainHess_claude.cpp`, `vm_file_structure.md` and the bookkeeping files, staged by path. No push without asking.

## Verification
1. The build links with no errors.
2. In the smoke test, Smd·conv equals the earlier `G_code_d` for the T2ss tissue to about 1e-9.
3. All of about 1007 `_SM.txt` files exist, and none are NaN.
4. The T2ss and ExTs Hessian means match the previous `G_vs_shapeindex_{T2ss,ExTs}.tsv` values.
5. The benchmark means printed by the script match `p0sm*` exactly.
