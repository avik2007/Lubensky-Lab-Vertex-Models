**Status:** APPROVED 2026-10-22 (supersedes the stub `2026-10-21_negative-mode-push.md`).

# Plan: negative-mode probe (detect → push → re-equilibrate → re-probe)

## Context
A Hessian eigenvalue λ < −threshold means the "equilibrated" tissue sits at a saddle, so its Hessian G is not a shear modulus.
The library's G function (`calculateShearModulusWithEigenDecompositiond`, `AGPerimVertex/AGPerimTissue.cpp:3005`) never logs these
modes: they go to the floppy test and are dropped or return 0. The main customer is QZ20 (QUARTIC2.0 double well), but the probe is
independent of QUARTIC2.0. User decisions (2026-10-22):
- A dedicated main file, used for all double-well Hessian runs, writes the negative-mode log.
- Zero mode: |λ| < 1e-8 · max(|λ_min|, |λ_max|). Negative mode: λ < −that threshold. The threshold lives in the library.
- Push: Δr = ± MINLENGTH · u, where u is the most negative mode's unit eigenvector (Eigen-normalised, norm 1 over all 2·N_v
  coordinates; MINLENGTH = 0.001 from the input file). Try +u and −u and keep the lower energy. No retries within a round (the doubled amplitude is only used in round 2, below).
- New functions go in AGPerimTissue (or AGCell). Double precision only. T1s stay on during re-equilibration.
- Recovered tissues carry a flag k = number of pushes. No separate detection-only pass: build it all, then smoke-test.
- **Rounds (user, 2026-10-22): at most 2.** Each round tries +u and −u once, nothing else. Round 1 uses amplitude s₁ = MINLENGTH.
  Round 2 doubles it: s₂ = 2·MINLENGTH, along the most negative mode of the CURRENT state. A tissue goes to round 2 if, after round 1,
  either neither sign lowered E (nothing was pushed), or the re-equilibrated tissue still has a negative mode. Still negative after
  round 2 (or neither sign lowers E in round 2) → STUCK, recorded.
- How STUCK tissues enter G averages is deliberately left undecided. The ledger records everything (G, λ_min, outcome per round),
  and we decide after seeing real data (smoke test, then QZ20).

## Library changes (`AGVertexModel/AGPerimVertex/AGPerimTissue.h/.cpp`)
1. Move `REL_EIGENVALUE_THRESHOLD = 1e-8` out of the G function into one class-level constant used by both the G function and the
   new detector. The value is unchanged: max|λ| = max(|λ_min|, |λ_max|), so G results don't change.
2. `int FindNegativeModesd(ofstream* errorfile, double* lambdaMin, double* lambdaMax, double* threshold, VectorXd* uMostNegative)`
   - Builds `calculateHessiand` and runs `SelfAdjointEigenSolver` (eigenvalues sorted, so the most negative is column 0).
   - Returns the count of λ < −threshold and gives back the most negative mode's eigenvector.
   - Index map: vertex i (in `vertices` list order) ↔ coordinates (2i, 2i+1), the same mapping as `calculateHessiand`.
3. `double PushAlongModed(const VectorXd& u, double amplitude)`
   - For each vertex: `rnew = r + amplitude·(u_2i, u_2i+1)`, then `HelperEvolve_UpdatePosition()`, which wraps positions back into
     the box and updates the edge flags.
   - Then refreshes every cell's `SavedInfoForEvolve_CWVerts` and returns `CalculateFreeEnergy()`.
4. `int PushAlongMostNegativeModed(const VectorXd& u, double amplitude, double* E0, double* Eplus, double* Eminus)`
   - Pushes +s, then −2s (net −s), and keeps the lower energy (push +2s back if +u was lower).
   - If neither is below E0: undo the push and return 0, which means STUCK. Otherwise return ±1.
All new code gets verbose comments; nothing is written to file inside the library except via `errorfile`.

## New probe driver `AGPerimHessian_test/AGPerimHessian_test/MainNegMode_claude.cpp`
Copy of `MainHess_claude.cpp`: same globals, `ReadInputData`, CWVerts set-up, and no Evolve.
- Usage: `./agperim_negmode_claude <EQ file> <out folder> <round k>`. The amplitude is 2^(k−1)·MINLENGTH (k = 1, 2), and k = 3 means
  probe only (the final check after round 2).
- Each call is one pass:
  1. Set CWVerts and call `FindNegativeModesd`.
  2. Compute Smd with `calculateShearModulusWithEigenDecompositiond` and write `_SM.txt`, as MainHess_claude does.
  3. If there is a negative mode: call `PushAlongMostNegativeModed`, then write the pushed state to `<name>_push<k>.txt` with
     `WriteOutputDataAppend` (a complete input file).
  4. Append one row to `<name>_negmode.tsv`:
     `pass  input  n_neg  lambda_min  lambda_max  threshold  E0  E_plus  E_minus  sign  amplitude  Smd  outcome`
     `outcome` is `MIN` (no negative mode), `PUSHED`, `NOPUSH` (neither sign lowers E; round 1 → go to round 2 on the same state),
     or `STUCK` (still negative at the final probe, or NOPUSH in round 2). The recovery flag k is the round after which the probe
     first said MIN.

## Orchestration `batch_negmode.sbat` (new work folder `/scratch/.../claude_VM_runs/NegMode/`)
- One array task per EQ file. The task loops:
  1. Probe pass k.
  2. If PUSHED: re-equilibrate `<name>_push<k>.txt` with the existing chunked loop and `agperim_equilibration_claude_chunk`
     (same CHUNK/CAP/USR1 rules as `batch_EQ_chunked.sbat`).
  3. Probe its `_EQ.txt` as pass k+1.
- Stop at MIN, STUCK, or CAP without EQ (recorded as `NE_CAP`). Never go beyond round 2.
- Makefile target `agperim_negmode_claude`, modeled on `Quartic_C4fix/Makefile`, with the same g++ flags and library sources compiled
  in place.

## Verification (smoke test)
1. Regression: on 3 T2ss + 3 QuOw EQ tissues, the probe's Smd equals `agperim_hess_claude`'s Smd bit-for-bit (the threshold refactor
   changes nothing).
2. Push mechanics (independent of whether any negative mode exists): for the lowest few positive modes of one tissue, check
   E(r ± s u) − E0 against ½ λ s² with s = 1e-3. A match to a few digits confirms the vertex↔index mapping, the sign and the
   wrapping. Include a tissue with vertices near the box edge.
3. Detection: probe ~20 existing EQ tissues (Quartic_C4fix MTQ3 C4 = 10, fluid p0 ≥ 3.84, T2ss solids) and report n_neg and
   λ_min. If one has a negative mode, run the full push → re-equilibrate → probe loop on it.
4. STOP: hand the library + driver diff to the user for review before any production run or commit (user code review is mandatory).
5. Bookkeeping: todo, recentactions, `vm_file_structure.md` (new driver + sbat), report §, plan file in `vm_claude_plans/`.

## Revision 2026-10-24 (code review #1 with the user): translations, push size, exact reset

These changes replace the fixed push size and the eigenvalue cut of the original plan. The code comments in
`MainNegMode_claude.cpp` (header, "PUSH SIZE") and `AGPerimTissue.cpp` (`FindNegativeModesd`,
`PushAlongMostNegativeModed`) repeat this explanation next to the code.

**1. Translation modes are removed explicitly.** In a periodic box of fixed size, moving every vertex by the same
vector costs no energy. So the Hessian always has two exact zero modes, the rigid translations
t_x = (1,0,1,0,…)/√N_v and t_y = (0,1,0,1,…)/√N_v (entries 2i and 2i+1 are vertex i's x and y). Their computed
eigenvalues are round-off, |λ| ≈ 1e-14. `FindNegativeModesd` sorts the eigenvalues by |λ| and takes the two
smallest as the translations. It then checks each one with a dot product: the squared overlap with the
translation plane, (u·t_x)² + (u·t_y)², must be above 0.99 (exactly 1 for a pure translation). If either check
fails, the probe reports TRANSFAIL and does not push.
- Measured on the 18 smoke tissues + T2ss rng150 3.88: always passed (overlap 1). The third smallest |λ| was ≥ 1.5e-11, at least ~800× above the translations.

**2. Negative modes = every non-translation mode with λ < 0.** There is no eigenvalue cut any more. The shear-modulus
zero-mode cut (1e-8·max|λ|) is only logged. Whether a weak mode can be acted on is decided by the push size (3).

**3. Push size from a target energy drop.** Variables:
- u = unit eigenvector of the most negative mode. Vertex i moves by s·(u_2i, u_2i+1).
- λ = its eigenvalue, the curvature of E along u [energy/length²].
- s = push size, the length of the whole displacement [length].
- E = total free energy before the push [energy].
- ε = DBL_EPSILON = 2.2e-16. Energies of size |E| are resolved only to ≈ ε|E|.
- ΔE = the energy drop aimed for [energy].
- ⟨ℓ⟩ = mean edge length [length].

To second order, E(r + s·u) − E(r) = ½λs². Setting the drop equal to ΔE gives

  **s = √(2ΔE/|λ|)**, with **ΔE = 100·ε·|E|·4^(k−1)** (k = round),

which has units of length. In round 1 the drop is 100× the energy round-off, so the choice between +u and −u
is a real comparison and not noise. Round 2 makes ΔE 4× larger, which doubles s.

**Cap:** s ≤ 0.1·⟨ℓ⟩, so no push can cross an edge or cause a T1. If the s needed is larger than the cap, the
outcome is MARGINAL and the tissue is not pushed. Examples (round 1):
- E = 3e5, λ = −4e-7 → ΔE = 6.7e-9, s = 0.18.
- E = 5, λ = −4e-7 → ΔE = 1.1e-13, s = 7.4e-4.

**4. Exact reset.** Trial pushes are undone by restoring the saved start state the way Evolve's backtracker does
(`Update_Old_Mechanical_Variables` + the old edge flags + fresh clockwise lists), not by an opposite push.
`currentFreeEnergy`/`pastFreeEnergy` are set to the state that is kept.
- Bit-exact on the tests.
- The user closed the boundary-wrapping test: the backtracker already handles wrapping.

**Test of the revised probe, round 1 (2026-10-24, scratchpad, 19 tissues):**
- Smd is bit-identical.
- Translations OK everywhere.
- Solids (p0 3.75 and 3.82): all MIN.
- Fluids (p0 3.86 to 3.88): 8 to 24 negative modes each (the old cut counted 4 of 18).
- Push outcomes:
  - PUSHED ×3: Aq2w and QuTs 3.86 (E ≈ 4), T2ss 3.88.
  - NOPUSH ×4: the quartic sets with E ≈ 1e5, where both signs raised E by 4–7× ΔE.

**Curvature scan along u**, using (E(r+su) + E(r−su) − 2E(r))/s². This is the energy's actual curvature at push size s:
- **The Hessian is correct.** The measured curvature → λ as s → 0 (T2ss 3.88: −3.43e-4 at s = 1e-3, λ = −3.429e-4).
- **The negative curvature only holds for small s.** At larger s the energy stiffens: curvature ≈ λ + c·s², with c ≈ 0.13 for the quartic fluids. So the energy decreases along u only up to s* ≈ √(|λ|/c):
  - Aq2w 3.86: s* ≈ 2e-3.
  - Lq2w 3.86: s* ≈ 8e-3.
- **The four NOPUSH tissues have s_needed > s*.** For the E ≈ 1e5 tissues, s_needed (0.018 to 0.03) is above s*, so the push lands where the energy is already rising. Their dip along u is only ~|λ|s*²/4 deep: about 1e-13 for Aq2w, about 1e-10 for Lq2w.
- **T2ss 3.88 is not at a stationary point.** It has a residual force along u of F·u ≈ 2e-3, so the push direction is set by that force.

**MANDATORY CODE REVIEW of this revision (user, 2026-10-24).** The revised push code (translation removal,
push size, exact reset, new outcomes MARGINAL/TRANSFAIL, the batch-loop column change) must be reviewed by the user
as a diff before it is committed or used for any run. Patch: `claude_VM_runs/Quartic2_Bidisperse/review/review1_negmode.patch`.

**Open (for the user):** how to treat NOPUSH. Doubling s in round 2 can only go further into the stiffening regime.
