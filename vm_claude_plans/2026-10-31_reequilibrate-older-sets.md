# Re-equilibrate the force-flagged older-set tissues

**Status:** IN PROGRESS. Approved by the user on 2026-10-31 ("All your guesses as far as decisions are good").
Job 65151016 (506 tasks) was submitted on 2026-10-31.

**Parent record:** [`2026-10-31_rk4-fallback-slow-tissue-incidents.md`](2026-10-31_rk4-fallback-slow-tissue-incidents.md).
That record covers the mechanism, the false equilibria and the residual-force check that found these tissues. The relevant
sections there are "Older sets: residual-force check" and "OPEN LEDGER". This plan carries out the user's decision recorded
in that file: re-scan all older sets and re-equilibrate whatever needs it.

## Background (summary; details in the parent record)
- The force check (job 65144451) covered 3139 older `_EQ.txt` files. 534 have a largest residual force Fmax at least 10× the
  median of their (set, p0) peers. Tiers: A 115 (≥ 100× and Fmax ≥ 1e-3), B 267 (≥ 30×), C 152 (10–30×).
- There are two causes:
  - 126 files were written at the time cap (T = 500/1000). The old drivers wrote `_EQ.txt` even when the run ran out of time.
  - 408 were declared EQ mid-run with large leftover forces. That is the false-EQ signature, consistent with the permanent fixed-step RK4 fallback.
- These EQ files feed the relthr Hessian G results used so far.

## Decisions (user, 2026-10-31)
1. **Scope: all flagged tissues, every tier,** not only tier A. Re-checking the forces afterwards also documents which tier-C flags were real.
2. **Output goes to scratch only:** `claude_VM_runs/ReEQ_older/`. The originals in `Avik/PhaseTransitionProject/`,
   `CompressionInducedFluidizationData/…` and `Quartic_C4fix/` are only read. Whether to replace them is the user's decision, later.
3. **QuOw / AqOw / LqOw (28 flagged) are skipped.** These are pre-C4-fix quartic sets, superseded by Qu2w/Aq2w/Lq2w.
   They are recorded in the ledger as `skipped_preC4fix`.
4. **Recompute G** for the re-equilibrated tissues and report the change in the per-(set, p0) averages.

## Method
- **Build:** the WARNING build (`agperim_equilibration_claude_chunk_warn`, which matches commit fe741fdf3), copied into `ReEQ_older/`.
- **Driver:** `batch_EQ_chunked_warn.sbat` with CHUNK=10, CAP=20000. Short chunks matter because every chunk is a new process,
  which resets the integrator to adaptive and the step to 1e-3. CAP is the cumulative time. The originals stopped at
  T ≤ 3946, so CAP is not a constraint.
- **Start state:** `make_reEQ_older.py` copies each flagged `<name>_EQ.txt` to `start/<set>/<name>.txt`. The driver keeps the
  file's CURRENTTIME, so times are cumulative. Tasks are ordered tier A, then B, then C, worst peer ratio first.
- **Smoke test** (2026-10-31, `ReEQ_older/smoketest/`): one file per set (25 sets), 0.02 time units each. Results:
  - All 25 parse and run with exit code 0, including the old formats with no C4/TISSUE_EXT_STRESS lines (the defaults are used).
  - Five flagged tissues (ExTs 100/3.79, L2Ex 120/3.77, L2Lo 100/3.75, LoTs 100/3.78, M2Lo 100/3.75) re-reached EQ within 0.014 time units.
  - Their Fmax fell from 0.02–0.31 to 3e-7–1.6e-5. The original large forces are real: recomputed and stored forces agree to 3e-11.
  - So the large leftover forces were localized and relax in about 0.01 time units of genuine dynamics. Those EQ verdicts were premature, not borderline.
- **Validation of the procedure on the Z6 false EQ** (job 65137664): the 18 seed-150 false-EQ tissues re-reached EQ.
  - Fmax fell from 1.1–7.8e-4 to 0.55–1.5e-5, the level of clean peers. Fmax/scale fell from about 6e-6 to about 1e-7.
  - Their final step sizes are about 1e-11 with `WARNING_RK4_FALLBACK`. At a true minimum, round-off in the energy rejects steps, so the WARNING alone does not mean false EQ. The force check decides.

## Steps
- [x] Build the candidate list (`forces_older/false_eq_candidates_tiered.tsv`).
- [x] Set up `ReEQ_older/` (`make_reEQ_older.py`): 506 queued, 28 skipped. Write `reEQ_older_ledger.tsv`.
- [x] Smoke test: all 25 sets.
- [x] Submit job 65151016 (array 1–506).
- [ ] When tasks finish: run the force check on every new `_EQ.txt` (`vmcheck_all.py`). Fill `status`/`Fmax_new` in the ledger.
      Compare against the peer medians.
- [ ] Tissues still NE at the job time limit: resubmit the same lines. They continue from `_NE.txt`.
- [ ] Recompute G (relthr Hessian, same settings as the original G) for the re-equilibrated tissues. Report the per-(set, p0)
      change in mean G, and the drop and add counts.
- [ ] Update the parent record (Actions, OPEN LEDGER) and this plan's Results section.
- [ ] Back up `ReEQ_older/` ledgers, the scripts and the new EQ files to turbo at session end. Scratch is purged.

## Files (`/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/ReEQ_older/`)
- `make_reEQ_older.py`: the setup script (copies only).
- `start/<set>/`: copies of the flagged original EQ files, used as starting states.
- `results/<set>/`: driver output (`_EQ.txt`, `_NE.txt`, `_error.txt`, `_chunks.tsv` with the rk4_fallback columns).
- `EQparameters_reEQ_older.txt`: 506 lines, `<start copy> results/<set>`.
- `reEQ_older_ledger.tsv`: one row per flagged file (534), with tier, original Fmax, peer ratio, cap flag, original path, status and new Fmax.
- `batch_EQ_chunked_warn.sbat`, `agperim_equilibration_claude_chunk_warn`: copies from `Quartic2_Sweep/`.
- `smoketest/`: the smoke-test script and its outputs.

## Results
(to be filled in)
