# Plan: fix the C4 shadowing bug, re-equilibrate the C4 ≠ 1 quartic sets, run the Hessian G on all quartic sets (no squeeze)

## Context
Before running the quartic sets, I found that forces, pressures and the Hessian have always used **C4 = 1**, while the energy used the file's C4.
- **Cause:** `AGTissue.h:41` declares a member `double C4;`. That member shadows the global `extern double C4` from `AGConstants.h`.
  - `ReadInputData` and `CalculateFreeEnergy` are tissue member functions, so they use the **member**, which holds the file value.
  - `AGCell::GetPressure` / `GetNegativeDPressureDArea` use the **global**, which every `Main*.cpp` leaves at 1.
- **Evidence:**
  - Every quartic file's pressure column implies C4 = 1.000.
  - QuTw, Qu2w, QuOw and MTQ3 have identical geometry. So do QuTs, Aq2w and AqOw, and so do Lq2w and LqOw.
- **Affected** (C4 ≠ 1): **QuTs (30), QuTw (3), Qu2w (2), Aq2w (2), Lq2w (2), MTQ3 (10)**. Their EQ tissues are really C4 = 1 tissues.
- **Unaffected:** QuOw, AqOw and LqOw (C4 = 1), plus all non-quartic sets.
- **User decisions:** fix the bug and redo the C4 ≠ 1 sets with the Hessian method only. **No squeeze runs.** A finite-difference sanity check is welcome.

Already checked analytically:
- With a = A/A0, E_A = ½KA0[(a−1)² + C4(a−1)⁴], Π = K[(1−a) + 2C4(1−a)³] and d²E/dA² = (K/A0)[1 + 6C4(a−1)²] are mutually consistent.
- The Hessian/G path calls the per-cell `GetPressure` and `GetNegativeDPressureDArea`, so the quartic potential works once C4 is passed through correctly.

## Step 0: the fix (repo, its own commit)
- In `AG Vertex Model/AGVertexModel/AGVertexModel/AGTissue.h`, delete the member `double C4;` (line 41) and leave a comment explaining why.
- Member functions then use the global C4, so `ReadInputData` sets the value that `AGCell` reads.
- No code uses `tissue.C4` explicitly (checked with grep).
- Rebuild `agperim_hess_claude`.

## Step 1: sanity checks (scratch: `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Quartic_C4fix/`)
**1a. Finite-difference check with the quartic energy**
- Copy the earlier FD driver (`T2ss_SM_Hess_testing/src/MainHessFDCheck.cpp` + `RefModel.h`) and add the quartic term to the reference energy:
  - E_A += ½(K/A0)·C4·(A−A0)⁴/A0²
  - the gradient weight wA gains + 2K·C4·(A−A0)³/A0³
- Run it on a few QuOw tissues and on one or two of the new C4 = 3 / C4 = 30 tissues from step 2's pilot.
- **Pass criteria:**
  - the energy matches `CalculateFreeEnergy`;
  - the library Hessian matches the FD Hessian (~1e-10);
  - Ξ matches FD;
  - G_code = G_FD;
  - forces are ~0 at the new equilibria.
- **Also confirm the fix:** the library pressure on a C4 = 3 file equals K[(1−a) + 6(1−a)³].

**1b. Equilibration regression (cheap)**
- On ~7 QuOw tissues (C4 = 1, where the bug changes nothing), check that the repo code re-equilibrates from the raw `_ELS_3.txt` to the existing `_EQ.txt` (same finalEnergy and positions).
- This confirms the repo library, which has diverged from your 2023 Great Lakes copy, still evolves the same way.
- **Gate:** if it doesn't, stop and report.
- Also check that the T2ss smoke tissue gives the same Smd after the fix (non-quartic results unchanged).

## Step 2: production (SLURM, `--account=abmondal1`, standard partition, arrays as before)
**Driver.** `AGPerimVertex/MainAGPerim_claude.cpp`: a minimal-change copy of repo `MainAGPerim.cpp` (currently `cin` + a Windows path) with argv input and an output folder.
- It follows your Great Lakes `MainAGPerim.cpp`: evolve until equilibrated or TOTALTIME, then write `_EQ.txt` with the initial/final energy footer.
- Makefile modeled on `GreatLakes_makefiles/makefile`.

**2a. Equilibrate the 6 C4 ≠ 1 sets from their raw `_ELS_3.txt`** (112 each, read-only):
- `PhaseTransitionProject/{QuTs,QuTw,Qu2w}_files`
- `CIF_HomeFolderBackupCouldHaveDuplicates/{Aq2w,Lq2w,MTQ3}_files`
- Pilot a few tissues first to get timing (your originals requested 140 h); report it before the full arrays.

**2a'. Overwrite the old EQ files (the user asked for this)**
- The equilibration runs write to scratch first. Only after the step 1 checks pass do the new `_EQ.txt` files overwrite the old C4-bug ones, in the folders where the old ones live:
  - `PhaseTransitionProject/{QuTs,QuTw,Qu2w}_files/`
  - `CIF_HomeFolderBackupCouldHaveDuplicates/{Aq2w_smaller,Lq2w_smaller,MTQ3}_files/`. The Aq2w and Lq2w EQ files live in `_smaller_files`; their raw inputs are in `_files`.
- **Safety copy first:** all old `_ELS_3_EQ.txt` files of the 6 sets (~47 MB) go to `/nfs/turbo/.../claude_hessian_testing/Quartic_C4fix/original_EQ_C4bug/`, keeping the folder structure.
- **Tissues that fail to equilibrate:** their stale C4-bug EQ file is moved into that safety folder, not left in place, so no bugged EQ remains. They are listed in the report.
- **Not touched:**
  - the old squeeze files (`_EQ_0.000N.txt`, `_EQ_0.000N_EQ.txt`);
  - the `_dE` folders;
  - duplicate copies of these sets in other folders (`CIF_VertexModelSimulationInputOutpudEFiles`, and the CIF backup's QuTs/QuTw/Qu2w copies). These are listed in the report.

**2b. Hessian G** (`MainHess_claude`, rebuilt with the fix, no evolution):
- the new EQ tissues of the 6 C4 ≠ 1 sets;
- the existing EQ files of QuOw, AqOw and LqOw, which are valid as they are.

## Step 3: plots and report
- **Figures:** your quartic overlays, rebuilt with the Hessian G and your markers/colors:
  - `T2ss_QuTw_Qu2w_QuOw` (C4 = 0, 1, 2, 3)
  - `Qu2w_QuOw`
  - `QuTw_T2ss`
  - `QuOw_AqOw_LqOw`
  - `Qu2w_Aq2w_Lq2w`
  - MTQ3 & QuTs
  - T2ss uses the existing Hessian results.
- **Squeeze comparison:** only for the C4 = 1 sets (QuOw/AqOw/LqOw), whose existing `.mat` files are valid.
- **Report:**
  - G per set;
  - equilibration success counts;
  - whether C4 now changes the geometry and G as expected.
- **Housekeeping:**
  - no original data is modified;
  - turbo backup (without the `_error.txt` files);
  - todo and recent-actions updates;
  - add the new driver to `vm_file_structure.md`.

## Commits (no push without asking)
1. `AGTissue.h` fix alone.
2. `MainAGPerim_claude.cpp`, `vm_file_structure.md` and the bookkeeping files.

## Verification summary
- 1a FD check passes on quartic tissues, including C4 = 3 and C4 = 30.
- 1b: the repo code reproduces your QuOw equilibria.
- T2ss Smd is unchanged after the fix.
- New C4 ≠ 1 tissues have pressure columns consistent with their header C4.
- QuTw, Qu2w and QuOw are no longer geometrically identical.
