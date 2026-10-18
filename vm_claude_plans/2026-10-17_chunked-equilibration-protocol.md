# Protocol: chunked, restartable equilibration (2026-10-17)

**Status (2026-10-17): APPROVED by the user. Cap = 10000. C1: keep running and restart if they don't finish. The user prefers "no restart" for C2 if it is possible.**
- Done: the driver (`agperim_equilibration_claude_chunk`) and `batch_EQ_chunked.sbat` are built. All §4 tests passed, locally and in a 10-min SLURM test (job 64217352: USR1 arrived, state saved at t = 7.618, job COMPLETED).
- Group A (3 MTQ3) is submitted as job **64217475**, starting from `_NE.txt` at t = 200 (originals kept as `_NE_leg1_t100.txt` and `_NE_NE.txt`).
- C2 was NOT cancelled. With the user's approval, 36 tasks were stopped with gdb (`set var tissue.equilibrated = 1`), which made the old driver write its state as `_EQ.txt`; `gdbstop_to_NE.sh` turned each into `_NE.txt` (original kept as `_gdbstop_EQ.txt`). They continue as job **64218833**. The 3 C1 tasks still run in 64105735.
Work folder: `/scratch/dkluben_root/dkluben1/abmondal/claude_VM_runs/Quartic_C4fix/`.

## Why
- The raw inputs set TOTALTIME per set: QuTs 1000, QuTw/Qu2w 10000, Aq2w/Lq2w/MTQ3 mostly 100000 (some 50000; 3 MTQ3 files at 100). These values come from the user's 2023/2025 files. The 2023 EQ files record them too, and those runs equilibrated at CURRENTTIME of about 3800–8900.
- SLURM kills a job at 72 h. With these TOTALTIMEs, a slow tissue times out before it reaches TOTALTIME. The current binary then writes nothing, so all of that compute is lost.
- The user's rule: keep each run short (TOTALTIME = 1000), always save the final state, and restart from it. Restarting is cheap. A run that produces no output is waste.

## 1. Driver changes (`AGPerimVertex/MainAGPerim_claude.cpp`)
1. **Chunk length from the command line.** `./agperim_equilibration_claude_NE <input> <outdir> [chunk]` runs for `chunk` time units (default: the file's TOTALTIME), overriding the file's value. The protocol always passes **1000**.
2. **Save on timeout.**
   - The driver catches SIGUSR1 and SIGTERM and sets a flag.
   - The main loop checks the flag after every `Evolve(TIMEINTERVAL)` call (one call is 0.001 time units, so the delay is short).
   - When the flag is set, it writes `_NE.txt` with the footer "stopped by signal", then exits.
   - SLURM sends USR1 15 min before the time limit (`--signal=B:USR1@900`). SIGTERM covers a cancel or the final kill.
3. **Cumulative time.**
   - CURRENTTIME is taken from the input file: 0 for a raw file, or the saved value for an `_NE.txt`.
   - The output records CURRENTTIME = start + time run in this chunk, so `_NE.txt`/`_EQ.txt` always show the tissue's total evolved time.
   - Each chunk still runs `chunk` time units from its own start. This part needs checking first: ReadInputData must read CURRENTTIME back, and nothing in Evolve may depend on its absolute value (mitosis is off).
4. **Stable file names.**
   - Outputs are always `<tissue>_ELS_3_EQ.txt` or `<tissue>_ELS_3_NE.txt`, no matter how many chunks ran. The driver strips trailing `_NE` from the input name, so there are no `_NE_NE_NE` chains.
   - Before overwriting `_NE.txt`, the previous state is kept as `_NE_prev.txt`, so a bad chunk can always be undone.
   - `_error.txt` is appended to, with a header line for each chunk.
5. **Chunk ledger.** One line is appended per chunk to `<tissue>_ELS_3_chunks.tsv` with these columns:

   | Column | Meaning |
   |---|---|
   | chunk # | chunk number |
   | start t | cumulative time at the start |
   | end t | cumulative time at the end |
   | initial E | energy at the start |
   | final E | energy at the end |
   | final slope | energy slope at the end |
   | outcome | EQ / NE / SIGNAL |
   | wall min | wall-clock minutes |

   It gives the whole history of a tissue at a glance.

Nothing else changes: the evolution loop and the equilibrium criteria stay as they are (|slope| < 1e-15; step size < 1e-14).

## 2. Batch file (`batch_EQ_chunked.sbat`)
- One array task per tissue: `--time=72:00:00` and `--signal=B:USR1@900`.
- **Loop:** while there is no `_EQ.txt`, run one chunk of 1000 from the latest state (`_NE.txt` if present, otherwise the raw input).
- **Stop when** the tissue equilibrates, the signal arrives (state is saved), or the cumulative cap is reached (§5).
- The driver runs in the background and the shell forwards USR1 to it (`trap` + `wait`).
- Nothing resubmits itself, and no watcher runs. Unfinished tissues are resubmitted by hand with the same command, and they continue from their `_NE.txt`.

## 3. Which tissues and what happens to each
| Group | Tissues | Current state | Action |
|---|---|---|---|
| A | 3 MTQ3 (rng150, p0 3.88/3.89/3.9) | done, `_NE_NE.txt`, cumulative t = 200 | Rename to `_NE.txt` (keeping the original), then start the protocol: chunks of 1000 |
| B | 33 QuTs in job 64198251 (second leg) | running with the `_NE` binary, TOTALTIME 1000; finishes in ~10 h with output | Leave running. When done, normalize names (`_NE_NE` → `_NE`, `_NE_EQ` → `_EQ`) and join the protocol if not equilibrated |
| C1 | 4 tissues in job 64105735 with \|slope\| < 1e-12 (Aq2w rng130 3.81 at 3.6e-15; Aq2w rng100 3.83; Aq2w rng120 3.83; Aq2w rng130 3.83) | t ≈ 3500–4000, about 37 h left before the limit | Keep running: these are the likeliest to equilibrate soon. If one times out, it reruns from raw under the protocol |
| C2 | the other 37 in job 64105735 (\|slope\| 1e-12 to 4e-9, including the 2 slow tissues at t ≈ 130–150) | t ≈ 3500–4000, likely to time out at t ≈ 7500 with no output | **Cancel now** and restart from raw under the protocol. Runs are bit-reproducible, so the state at t ≈ 3800 comes back after ~4 chunks; only the progress since then is redone |
| D | 1 QuTs already `_NE_EQ.txt` | equilibrated at cumulative t = 2000 | Rename to `_EQ.txt` (keeping the original) |

## 4. Verification before submitting
1. **Signal test:** a short run under `timeout -s USR1`. It must write `_NE.txt`, and a restart from it must start at exactly that energy.
2. **Cumulative time:** a chunk restarted from an `_NE.txt` at t = T must record T + 1000.
3. **Names:** a raw input followed by two restarts must give one `_NE.txt`, one `_NE_prev.txt`, and a 3-line ledger. On equilibration it must give `<tissue>_ELS_3_EQ.txt`.
4. **Reproducibility:** chunking must not change the trajectory within a chunk (log identical to a straight run, as checked on 2026-10-16). The chunk boundary resets the slope window, energyList, and step size to INITIALSTEPSIZE. That is the user's normal restart behavior and is accepted.
5. **In SLURM:** a 10-minute test job with `--signal=B:USR1@300` to confirm the signal reaches the driver.

## 5. Open questions for the user
1. **Cumulative cap:** after how much total time should a tissue that has not equilibrated be stopped and reported? Proposed: 10 chunks (t = 10000). The 2023 runs needed up to ~8900.
2. **Group C1:** keep those 4 running, or cancel everything for a clean start?

## 6. Bookkeeping
- Update the todo, the report (§8 D2 and §11), and the README in the work folder.
- Commit the driver and the bookkeeping, staged by explicit path.
- The downstream steps (Hessian, overwriting the EQ files, plots) are unchanged; they use `<tissue>_ELS_3_EQ.txt`.
