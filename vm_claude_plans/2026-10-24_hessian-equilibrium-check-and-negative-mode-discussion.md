# Hessian check of equilibrium, and the negative-mode discussion (record for the paper)

**Status:** FILED AWAY — potential plan, NOT approved, do not implement until the user asks. Part 1 is the proposed
plan; Part 2 is the record of the 2026-10-24 discussion between Avik and Claude that led to it, kept so it can be
recalled when writing the paper.

Related: `2026-10-22_negative-mode-probe.md` (the probe, with its "Revision 2026-10-24" section),
`2026-10-17_chunked-equilibration-protocol.md` (the current equilibration protocol).

---

## Part 1 — Potential plan: add a Hessian check to the equilibration

**Idea (user, 2026-10-24).** The equilibration currently stops when the slope of the normalized energy vs. time,
dE/dt, falls below `equilibriumtolerance` (1e-13 in the drivers, 1e-15 in the PhaseTransitionProject files). It
never looks at forces or curvature. Keep that test exactly as it is, and add a second test after it passes:

1. dE/dt test passes (unchanged).
2. Build the Hessian of the final state and diagonalize it (`calculateHessiand`).
3. Remove the two rigid translations: the two modes with the smallest |λ|, each checked by its dot product with
   t_x, t_y (overlap > 0.99), as in `FindNegativeModesd`.
4. Look for negative eigenmodes that are not translations (λ < 0).
   - None → the state is a true minimum: write `_EQ.txt` as now.
   - Some → the state is not a minimum: record it (count, λ_min, residual force), and decide what to do
     (to be planned: continue relaxing / push along the mode / flag).

**MANDATORY CODE REVIEW (user, 2026-10-24).** Any code written for this check, and any push code it uses
(`FindNegativeModesd`, `PushAlongModed`, `PushAlongMostNegativeModed`, the push-size rule in `MainNegMode_claude.cpp`),
must be handed to the user as a diff and reviewed BEFORE it is committed or used for any run. That includes the push
changes made on 2026-10-24 (translation removal, s = √(2ΔE/|λ|) with cap, exact reset); they have NOT been reviewed yet.

**Not decided yet (to settle when this is picked up):**
- What happens to a tissue that fails step 4 (continue the evolution, push, or only flag it).
- Whether to also record the residual force (max |F| per vertex and F·u along the negative mode), which separates
  "not relaxed" from "relaxed but flat" (see Part 2, §5).
- The cost: one Hessian diagonalization (2N_v × 2N_v, ~400 × 400 here) per equilibrated tissue, small next to the
  equilibration itself.

---

## Part 2 — Record of the discussion (2026-10-24)

### 1. Separating translations from other modes
In a periodic box of fixed size, moving every vertex by the same vector costs no energy, so the Hessian always has
two exact zero modes, the rigid translations

  t_x = (1,0,1,0,…)/√N_v,  t_y = (0,1,0,1,…)/√N_v   (entries 2i, 2i+1 = vertex i's x, y).

Numerically their eigenvalues are round-off, |λ| ≈ 1e-14 (|H·t| ≈ 3e-14). User's proposal, adopted: sort the
eigenvalues by |λ|; the two smallest are the translations; confirm each with a dot product, (u·t_x)² + (u·t_y)² > 0.99.
Measured on 19 tissues (solid and fluid, quadratic and quartic): always passed, overlap = 1; the third smallest |λ|
was ≥ 1.5e-11, i.e. at least ~800× above the translations.

### 2. Why not the determinant
Considered (user) and dropped (agreed): det H = product of all eigenvalues. It is always 0 because of the two
translations (numerically its sign is the sign of round-off); even with them removed its sign only gives the parity
of the number of negative modes; and a product of ~400 numbers between ~1e-5 and ~600 overflows/underflows a double.
The eigenvalue list already contains strictly more information.

### 3. Push size along a negative mode (dimensionally consistent form)
Variables: u = unit eigenvector of the most negative mode (vertex i moves by s·(u_2i, u_2i+1)); λ = its eigenvalue,
the curvature of E along u [energy/length²]; s = push size, the length of the whole displacement [length];
E = total free energy [energy]; ε = DBL_EPSILON = 2.2e-16 (energies of size |E| are resolved only to ≈ ε|E|);
ΔE = target energy drop [energy]; ⟨ℓ⟩ = mean edge length [length].

To second order E(r + s·u) − E(r) = ½λs², so a drop ΔE needs

  **s = √(2ΔE/|λ|)**  [√(energy / (energy/length²)) = length],  ΔE = 100·ε·|E| (round 1), capped at s ≤ 0.1⟨ℓ⟩.

(Earlier options written as bare scalings "1/|λ|", "1/√|λ|" were rejected by the user because their units are wrong;
e.g. MINLENGTH/|λ| is not a length.) The old fixed push s = MINLENGTH = 1e-6 gave drops of ~1e-16, below the energy
resolution, so the choice between +u and −u was round-off noise.

### 4. What the curvature scan showed
Finite-difference curvature along u at push size s: (E(r+su) + E(r−su) − 2E(r))/s².
- **The Hessian is right.** As s → 0 the curvature → λ (T2ss rng150 p0 3.88: −3.432e-4 at s = 1e-3 vs λ = −3.429e-4).
- **The energy stiffens quickly** in the fluid tissues: curvature ≈ λ + c·s², c ≈ 0.13 (quartic sets), so the energy
  only decreases along u up to s* ≈ √(|λ|/c) (≈ 2e-3 for Aq2w 3.86, ≈ 8e-3 for Lq2w 3.86). The dip along u is
  only ~|λ|s*²/4 deep (~1e-13 to 1e-10). For tissues with large E (~1e5, quartic sets) the push needed to resolve the
  energy (0.018–0.03) is larger than s*, so both signs raise E ("NOPUSH"). These look like floppy fluid modes held
  by higher-order terms, not real saddles.
- Round-1 probe on 19 tissues: solids (p0 3.75, 3.82) all minima; fluids (p0 3.86–3.88) have 8–24 negative modes each.

### 5. T2ss rng150 p0 3.88: a negative eigenvalue means the state is not a minimum
User's point: at a true minimum the lowest eigenvalues are the translations (≈ 0) and everything else is ≥ 0; a
non-translation eigenvalue below zero means equilibration did not fully succeed. Agreed. This tissue has 8 such modes,
λ_min = −3.4e-4 (10 orders above the translation round-off in magnitude), and a residual force along u of
F·u ≈ 2e-3 (constant over s = 1e-4 to 0.1): it is not even a stationary point. It passed the equilibrium test because
that test is on dE/dt, not on forces or curvature (for the quartic re-equilibration, 544/607 tissues passed with the
step size collapsed and the last energies identical; not checked for this T2ss file).

Two different situations therefore show up as "negative mode":
- **not relaxed** — non-negligible force and negative curvature (T2ss 3.88);
- **relaxed but flat** — force ≈ 0 (F·u ~ 1e-8) and a negative dip too shallow to act on (quartic fluids at p0 3.86).

**User decision:** not worried about T2ss — the Hessian shear modulus for T2ss agrees with the squeeze benchmark and
with the other potentials, so no residual-force study of the T2ss files. Instead, file the equilibration check of
Part 1 for later.

### 6. Still open from this discussion
- How the probe should treat NOPUSH (proposed: final outcome "SHALLOW" instead of round 2 with a doubled push, since
  a bigger push only goes further into the stiffening). Not decided.
- Code review #1 (negative-mode probe) and #2 (QUARTIC2.0 per-cell A0/P0 energy) are awaiting the user's approval.
