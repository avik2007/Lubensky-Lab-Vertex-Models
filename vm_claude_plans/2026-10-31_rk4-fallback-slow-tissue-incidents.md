**Status:** RECORD for the paper (not a plan to execute). Opened 2026-10-31. Every slow-tissue / throttled-step incident goes here,
with the date, the tissue, the evidence and what was done. The user asked for this explicitly ("I'll need them for writing a paper").

# Slow-tissue incidents: the fixed-step RK4 fallback in `AGPerimTissue::Evolve`

## Mechanism (diagnosed 2026-10-31)
- In `AGPerimTissue.cpp` (around line 691 of the committed library; the sweep binaries are built from `src_head/`), a step that
  RAISES the free energy is rolled back, the step size is multiplied by 0.75, and the global `INTEGRATORNAME` is set to `"RK4"`
  (fixed step). Nothing ever sets it back to `"RK4_Adaptive"`.
- So after the first energy increase, the step can only shrink, by 0.75 at each later increase, until the process exits.
  A new process (each chunk of `batch_EQ_chunked.sbat`) resets the integrator to adaptive and the step to INITIALSTEPSIZE = 1e-3.
- The largest possible step is 1e-3, because TIMEINTERVAL = 1e-3 clips the adaptive step.
- The box plays no part. TISSUE_ALPHA = 2.8e18 freezes it: the box error estimate is exactly 0 and LX = LY = 353.553 all along.

## Evidence
- AqZ6 rng 150 p0 3.81 (job 64834279_?, 72 h wall time): dt held at 7.617e-5 from t ≈ 4.3, then at 1.3557e-5 for the rest of the run.
  The ratio is exactly 0.75^6.000. After 4304 wall minutes it had reached only t = 268 and never finished chunk 1 (1000 time units).
  The `_error.txt` is about 5 GB. L2Z6 and M2Z6 at the same seed and p0 stalled the same way (t = 246 and t = 261).
- Debug restart of the AqZ6 one from its `_NE.txt` (2026-10-31, `Quartic2_Sweep/check/stepdbg/`, debug-only copy of the library
  that prints the error-limiting component): it ran at dt = 1e-3 from the first step. The largest error was about 1e-10 of the
  tolerance and came from a vertex; the box error was 0. So the state itself does not need small steps.
- AqZ6 rng 150 p0 3.80: declared EQ at t = 42.8 while throttled (dt levels 3.117e-5 then 4.161e-6, ratio exactly 0.75^7.000).
  **An EQ declared while throttled is suspect.** The dE/dt test is made over a few tiny steps. Check its forces before using its G.

## Detection (`Quartic2_Sweep/analysis/slow_steps/scan_steps.py`)
- Samples 300 windows per log, evenly spaced in bytes (so evenly in steps), and measures dt from consecutive step times.
- `fixed_small` = dt < 0.99e-3 and constant within a window. This is an inference: the adaptive controller can also hold dt
  constant. The support is that successive constant levels differ by exact powers of 0.75 (column `level_ratios_as_pow075`).
- `frac_steps_fixed_small` is the fraction of the tissue's steps (roughly its wall time) spent throttled. `throttled = yes` if ≥ 5 %.
- Outputs per sweep: `slow_step_ledger.tsv` (one row per tissue) and `slow_step_timeline.tsv` (dt versus t, for figures).

## Incident log
| date | set(s) | what | evidence / file | action |
|---|---|---|---|---|
| 2026-10-31 | L2Z6, AqZ6, M2Z6 rng 150 p0 3.81 | stalled for 72 h in chunk 1 (t ≈ 250–270) | see Evidence | held back from resubmission 65135510 |
| 2026-10-31 | AqZ6 rng 150 p0 3.80 | EQ declared while throttled (t = 42.8) | ledger row | flag; force check before G |
| 2026-10-31 | C4 = 1/6 sweep, all 336 | full scan (job 65136884): 114 throttled (53 EQ, 61 NE) | `analysis/slow_steps/Z6/` | see Scan results |
| 2026-10-31 | L2Z6, AqZ6, M2Z6 rng 150 p0 3.75–3.80 (18) | FALSE EQ: declared EQ at t ≤ 43 at dt ~1e-6; Fmax 10–60× clean | `Z6/forces_eq_throttled.tsv` | exclude from G; re-equilibration pending user decision |
| 2026-10-31 | 35 L2Z6/M2Z6 EQ tissues, mostly p0 ≥ 3.84 | mild throttle (dt ≈ 4e-4); Fmax about 2× clean | ledger rows | kept, flagged |
| 2026-10-31 | 61 NE tissues (58 AqZ6) | throttled while running out of wall time | ledger rows | resubmitted in 65135510 (step resets at each chunk) |

## Scan results
### C4 = 1/6 sweep, first pass (job 64834279; scan job 65136884, 2026-10-31)
Files: `Quartic2_Sweep/analysis/slow_steps/Z6/` (`slow_step_ledger.tsv`, `slow_step_timeline.tsv`, `forces_eq_throttled.tsv`,
`forces_eq_clean.tsv`). "Throttled" means at least 5 % of a tissue's steps were in the fixed-small-step mode.

| set | tissues | throttled | throttled & EQ | throttled & NE | median throttled fraction of steps (among throttled) |
|---|---|---|---|---|---|
| L2Z6 | 112 | 20 | 19 | 1 | 0.97 |
| AqZ6 | 112 | 64 | 6 | 58 | 0.64 |
| M2Z6 | 112 | 30 | 28 | 2 | 0.98 |
| all | 336 | 114 | 53 | 61 | |

- 61 of the 109 NE tissues (all stopped by the 72 h limit) were throttled; the other 48 were not.
  AqZ6 dominates: 58 of its 72 NE tissues were throttled.
- Step-level evidence: 108 jumps between constant dt levels. 67 are whole powers of 0.75 (34 down, 33 up). An upward jump is a
  reset at a chunk restart, followed by a new fall. The other 41 are restarts that settled on a level unrelated to the old one.
- **False equilibria: seed 150, p0 = 3.75–3.80, all three sets (18 tissues).** Declared EQ at t = 4–43 while dt ≈ 1e-6 to 5e-6
  (min 7.6e-7). Their residual Fmax is 1e-4 to 8e-4, compared with about 1.2e-5 for clean EQ tissues at the same p0 (for example rng 100 at 3.75).
  That is 10 to 60 times larger. Relative to the edge-tension scale Fmax/scale is about 5e-6, against about 1e-7 for clean tissues.
  These are NOT equilibrated. The dE/dt test passed only because the steps were tiny. Their G (job 65135521) must not be used.
  Seed 150 at p0 = 3.81 is the stalled NE triple above, so seed 150 is throttled at every p0 from 3.75 to 3.81.
- The other 35 throttled-EQ tissues (L2Z6/M2Z6, mostly p0 ≥ 3.84, seeds 100/110/120/130/160) had a mild throttle
  (dt ≈ 3.5–4.6e-4, about 0.75^3 of 1e-3; a few reached 1e-5). Their absolute Fmax is about 1.2–2.2e-5, compared with about 1e-5 for clean EQ tissues at the same p0.
  That is about 2 times larger, and it is not flagged as unequilibrated. Fmax/scale is large for them (up to 16) only because
  P ≈ P0 at high p0 makes the tension scale tiny (1e-6). Clean EQ tissues at those p0 show the same artefact (up to 15).
- AqZ6 throttled-EQ (6 tissues): Fmax/scale ≈ 6e-6, against about 1.5e-7 for clean AqZ6 EQ. Those 6 are seed 150, p0 3.75–3.80 (the false equilibria above).

### Older sets: residual-force check (job 65144451, 2026-10-31)
Files: `Quartic2_Sweep/analysis/slow_steps/forces_older/`
- `all_forces_with_peer_ratio.tsv`: every file.
- `false_eq_candidates.tsv`: tissues with Fmax at least 10× the median of their set and p0.
- `false_eq_candidates_tiered.tsv`: the same, with tier, cap and footer columns.

- 3139 EQ files checked, 0 parse errors. **534 flagged.**
- Tiers:
  - A (≥ 100× AND Fmax ≥ 1e-3): **115**. Robust: large in absolute terms too.
  - B (≥ 30×): 267.
  - C (10–30×): 152. At high p0 the peer median can be about 1e-6, which inflates the ratio, so tier C is the least certain.
- **126 of the 534 were written at a round time cap (T = 500 or 1000).** The old drivers wrote `_EQ.txt` even when TOTALTIME
  ran out (see the `MainAGPerim_claude.cpp` header). So these are unequilibrated states labelled EQ, a separate cause from the
  RK4 fallback. The other 408 were declared EQ mid-run with large leftover forces, which is the false-EQ signature.
- Flagged per set (flagged / checked): Aq2w 17/120, AqOw 11/106, ExTs_LP 37/112, ExTs 15/112, L2Ex_LP 30/110, L2Ex 10/112,
  L2Lo 11/112, L2Lo_LP 29/110, L2Ts 15/112, L2Ts_LP 37/111, LoTs 11/112, LoTs_LP 35/112, Lq2w 6/114, LqOw 8/111, M2Ex_LP 32/112,
  M2Ex 13/111, M2Lo 22/112, M2Lo_LP 33/112, M2Ts 15/112, M2Ts_LP 40/112, MTQ3 5/113, Qu2w 6/118, QuOw 9/107, QuTs 9/112,
  QuTw 6/116, SWAB 19/112, T2ss 18/112, T2ss_LP 35/112. The _LP (high-polydispersity) sets are hit about 2–3× as often.
- By seed: 150 → 244, 160 → 113, 130 → 92, 110 → 30, 100 → 29, 120 → 17, 140 → 9. Seed 150 dominates, as in the C4 = 1/6 sweep.
- Worst cases include T2ss_LP seed 150 p0 3.89 (Fmax 0.16, about 5.6e4× its peers), AqOw seed 160 p0 3.89, T2ss seed 160 p0 3.89
  (T_eq ≈ 28, 0.1× peers) and M2Ex_LP seed 160 p0 3.83 (Fmax 0.58).
- 51 files have stored forces that differ from the recomputed ones by more than 1e-6, which means stale force columns. This is a separate note.
- **These EQ files feed the Hessian G results (relthr) used so far.** Flagged tissues should be left out of G, or re-equilibrated.
- **Re-equilibration of the flagged older tissues:** plan [`2026-10-31_reequilibrate-older-sets.md`](2026-10-31_reequilibrate-older-sets.md),
  job 65151016 (506 tissues; QuOw/AqOw/LqOw skipped). Smoke test: in five tier A/B tissues, Fmax 0.02–0.31 relaxed to ≤ 1.6e-5
  within 0.014 time units, so those leftover forces were localized and the EQ verdicts premature.

### Quartic_C4fix sets: log scan (job 65137748) cross-checked with forces
| set | throttled (all EQ) | force-flagged | both |
|---|---|---|---|
| Aq2w | 20 | 17 | 7 |
| Lq2w | 22 | 6 | 6 |
| MTQ3 | 20 | 5 | 5 |
| Qu2w | 19 | 6 | 6 |
| QuTs | 21 | 9 | 7 |
| QuTw | 18 | 6 | 6 |
37 of the 49 force-flagged tissues were throttled. This is independent support that the RK4 fallback causes the false EQ.

## User decisions (2026-10-31)
- NO library fix: the fallback behaviour stays as it is.
- A WARNING must appear in the outputs whenever RK4_Adaptive switches to RK4, so individual users can decide whether to re-run.
- Re-equilibrate the 18 false EQ tissues and the 61 throttled NE tissues.
- Re-scan ALL older sets. Run as many as Great Lakes allows and keep an OPEN LEDGER, below, of tissues that need re-equilibration.

## WARNING output (added 2026-10-31; output only, the dynamics are unchanged)
- `AGPerimTissue::Evolve` writes to `_error.txt` at the first switch:
  `WARNING: RK4_Adaptive -> fixed-step RK4 at chunk time t = …, step size now …`. The same line goes to stdout, so it reaches the SLURM .out log.
  Every later shrink writes `WARNING: fixed-step RK4 step size reduced to … at chunk time t = …`.
- `MainAGPerim_claude.cpp` adds a footer line to `_EQ.txt` / `_NE.txt`:
  `WARNING: RK4_Adaptive switched to fixed-step RK4 during this run (final step size …). … consider re-running from this file.`
  It also adds two columns at the end of each `_chunks.tsv` row: `rk4_fallback` (`WARNING_RK4_FALLBACK` or `none`) and `final_stepsize`.
  Ledgers created before this change keep their old 8-column header.
- Repo: COMMITTED 2026-10-31 as fe741fdf3 (user OK). Only the warning hunks were staged; the uncommitted review-#1/#2 edits
  in `AGPerimTissue.cpp/.h` stay in the working tree. The committed code is byte-identical to the `src_warn` build that ran.
- Sweep build: `Quartic2_Sweep/src_warn/` (= `src_head/` + the patch), `Makefile.warn` → `agperim_equilibration_claude_chunk_warn`,
  run by `batch_EQ_chunked_warn.sbat`. The old binaries are untouched (job 65136026 is still using them).

- Verified 2026-10-31. In the re-equilibration runs (job 65137664), the seed-150 tissues print
  `WARNING: RK4_Adaptive -> fixed-step RK4 at chunk time t ≈ 1.6, step size now 0.00375`.
- Finding: the adaptive step VARIABLE can grow past 1e-3 (here to 0.005), even though TIMEINTERVAL clips each adaptive step to
  1e-3. A fallback from such a step leaves a fixed step LARGER than 1e-3 (0.75 × 0.005 = 0.00375). That only throttles after
  more energy increases. The scan counts dt ≥ 0.99e-3 as "full", so it reports throttling, not every fallback event.
  The WARNING lines record every event.

## Actions
| date | what | job |
|---|---|---|
| 2026-10-31 | 58 throttled NE tasks of resubmission 65135510 stopped with USR1. All 58 saved `_NE.txt` (outcome SIGNAL). | 65135510 |
| 2026-10-31 | 18 false EQ: `_EQ.txt` → `_EQ_false_2026-10-31.txt` (kept) and copied to `_NE.txt` to restart from that state. Their G moved to `relthr_falseEQ/<code>/`. | — |
| 2026-10-31 | Re-equilibration of 79 tissues (18 false EQ + 61 throttled NE, including the 3 seed-150 p0 3.81) on the WARNING build, CHUNK=10, CAP=20000. `EQparameters_Z6_reEQ_throttled.txt` | 65137664 |
| 2026-10-31 | Residual-force check of 3139 older EQ files (22 original sets + the 6 Quartic_C4fix sets). The user ran the originals and they have NO logs, so only false equilibria (leftover forces) and time-to-EQ can be checked, not throttling. Checker `analysis/slow_steps/vmcheck_all.py` (adds LOG and EXPONENTIAL; it matches the file forces to 2e-11 and passes the FD check on all 5 energy types). Output `analysis/slow_steps/forces_older/`. First hint: L2Lo_EQ rng 100 p0 3.75 has Fmax 1.8e-2, about 1000× the others. | 65144451 |
| 2026-10-31 | Scan of the 6 Quartic_C4fix sets (Aq2w, Lq2w, MTQ3, Qu2w, QuTs, QuTw; 672 logs). Output `analysis/slow_steps/<set>/` | 65137748 |
| 2026-10-31 | Partial result of 65137664: 19 of 79 have re-reached EQ (all 18 seed-150 false EQ, plus L2Z6 150/3.81). Fmax fell from 1.1–7.8e-4 (false EQ) to 0.55–1.5e-5, the level of clean peers. Their final steps are about 1e-11 with WARNING_RK4_FALLBACK, from round-off at a true minimum, so the WARNING alone does not mean false EQ; check the forces. | 65137664 |
| 2026-10-31 | Re-equilibration of the 506 force-flagged older tissues (all tiers; QuOw/AqOw/LqOw skipped) into `claude_VM_runs/ReEQ_older/`. Originals untouched. Plan `2026-10-31_reequilibrate-older-sets.md`. | 65151016 |

## OPEN LEDGER: sets and tissues needing re-equilibration
| set | logs at | scanned | throttled | false EQ (force check) | re-equilibration |
|---|---|---|---|---|---|
| L2Z6 / AqZ6 / M2Z6 | Quartic2_Sweep/results | yes | 20 / 64 / 30 | 6 / 6 / 6 | job 65137664 |
| L2Z1 / AqZ1 / M2Z1 | Quartic2_Sweep/results | no: sweep STOPPED by the user 2026-10-31 (150 EQ, 186 NE saved); scan the logs if it is ever used | | | not planned |
| Aq2w, Lq2w, MTQ3, Qu2w, QuTs, QuTw | Quartic_C4fix/results | job 65137748 | | | |
| 22 original sets: L2Ts/T2ss/M2Ts, L2Lo/LoTs/M2Lo, L2Ex/ExTs/M2Ex, QuOw, SWAB, the _LP versions (`Avik/PhaseTransitionProject/<set>_EQ` or `_files`); LqOw/AqOw (`CompressionInducedFluidizationData/CIF_HomeFolderBackupCouldHaveDuplicates/`) | no logs (user, 2026-10-31) | n/a | n/a | force check job 65144451: 534 flagged (incl. 6 C4fix sets) | job 65151016 (506; 28 *Ow skipped), ledger `ReEQ_older/reEQ_older_ledger.tsv` |

## Open questions
- How many older sets (Quartic_C4fix, Hess_SM_groups sources, L2Ts/T2ss/M2Ts …) were throttled the same way? They were
  equilibrated with drivers that share this `Evolve`. Their logs, where they still exist, can be scanned the same way.
- Library fix: DECLINED by the user on 2026-10-31. The workaround is short chunks plus the WARNING output.

## 2026-11-05 update: round 2 finished, round 3 submitted
- Round 2 (jobs 65135510, 65137664, 65151016) all ended. Most tasks were stopped by the 72 h signal, not equilibrated.
- Z6: 242 EQ / 94 NE (L2Z6 9, AqZ6 69, M2Z6 16 NE; all NE at p0 >= 3.79). 5 AqZ6 tissues reached CAP 20000 (`EQparameters_Z6_capped.txt`): capped, excluded from G.
- Older sets: 441 EQ / 65 NE (NE spread over 18 sets, t ~ 7100-11600).
- Round 3 (CHUNK=10, CAP=20000, WARNING build): job 65462965 (89 Z6), job 65462966 (65 older). By the user's rule, tissues still NE after round 3 are investigated, not resubmitted. Every one of them is to be recorded here.
