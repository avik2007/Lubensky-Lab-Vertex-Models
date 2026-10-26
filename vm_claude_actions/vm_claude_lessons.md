# Claude Lessons Learned

## Format
Each lesson: **Mistake** → **Lesson** → **How to Apply**

---

## Lessons

### L001 — Don't assume you know a named reference
**Mistake**: When the user mentioned "Boris Chernsky's recommended CLAUDE.md", I searched the web and the repo but couldn't find it, yet proceeded to create files without the template. I should have asked immediately.
**Lesson**: If a task depends on a named resource (person, document, template) you cannot locate, ask before doing anything that will need to be redone.
**How to Apply**: Before building on any named external reference, confirm you have it in hand. Don't guess or substitute.

### L002 — "Ignore language X" means ignore anything that only serves X, not just files of that extension
**Mistake**: A prior plan (2026-10-09) kept a C++ fix (`CalculateShearModulusForPythonCalculation`, which exports a Hessian purely so a Python script can eigendecompose it) in scope, reasoning "the Python script itself is excluded but this C++ function isn't the script." The user corrected this 2026-10-10: Python was only ever a quick side attempt, and now that Claude is doing the work the whole pipeline stays in C++ — so any code whose sole purpose is to feed the excluded path is excluded too.
**Lesson**: When a user puts a language/tool out of scope, check for same-repo code in an in-scope language whose only reason to exist is to serve that out-of-scope path, and exclude it as well. Don't scope by file extension alone.
**How to Apply**: Before including a function/file in scope, ask "does this exist for its own sake, or only to feed the excluded thing?" If the latter, leave it out unless the user asks for it specifically.

### L003 — At close-out, commit all of the session's bookkeeping together
**Mistake**: (2026-10-14) Closing out, I listed three uncommitted items: the todo, the recent-actions log, and the untracked `vm_claude_plans/`. The user said "commit that folder", and I committed only `vm_claude_plans/`. The user wanted the `vm_claude_actions` logs committed too and had to ask again, frustrated.
**Lesson**: The session's agent-state files (`vm_claude_actions/*`, `vm_claude_plans/*`) form one unit. When the user approves committing part of it at close-out, the rest of that bookkeeping belongs in the same pass.
**How to Apply**: When asked to commit at close-out, include every changed `vm_claude_actions/` and `vm_claude_plans/` file from the session, staged by explicit path, never `git add -A` (the CRLF drift). Ask only if source code or data would be included.

### L004 — Report the quantity the user asked for; don't push my own side-variants
**Mistake**: (2026-10-14/15) I kept attaching G_pure and ½(G_simple + G_pure) to every G comparison and framed the next step as a choice about them. The user said: "DON'T GIVE ME THIS NONSENSE ABOUT G_SIMPLE+G_PURE."
**Lesson**: The deliverable is the library's Hessian G (simple shear) against the squeeze benchmark. Derived variants I came up with during a diagnosis are not deliverables.
**How to Apply**: Plot and report only the series the user names. Mention a side-variant at most once, and only if it changes a conclusion. Use the x axis the user asks for (here, the shape index P0/sqrt(<A>)).

### L005 — "Fixed" means the user sees no gaps: check every figure and table for holes, and don't mask data silently
**Mistake**: (2026-10-15) After switching floppy tissues to G = 0, I reported the holes as fixed. But I had deliberately cut the ratio panel of `C4one_vs_T2ss` at p0 = 3.81 (0/0 in the fluid), and the user had stale downloads. The user: "DID YOU ACTUALLY FIX THE PROBLEM, I STILL SEE HOLES IN THE PLOTS". I had also first explained the floppy discrepancy as "absolute thresholds scale with K" before checking. The real cause was the LU-rank loop start plus round-off.
**Lesson**: A data-completeness fix is not done until every affected table has been checked for NaN/blank entries and every affected figure has been viewed. Any masking I add on purpose (cut ranges, dropped points) must be said up front, or better, replaced with a quantity that is defined everywhere (e.g. a difference instead of a ratio). Check a diagnosis against the data before stating it.
**How to Apply**: After any change meant to remove gaps, grep the summary tables for nan/blank, Read each regenerated PNG, and tell the user to re-download files that changed. Present a mechanism as a hypothesis until a measurement confirms it.

### L006 — Before writing a new driver, look for the user's own version of every behavior it needs
**Mistake**: (2026-10-15/16) I wrote `MainAGPerim_claude.cpp` from the user's `~/MainAGPerim.cpp` and added "write `_EQ.txt` only if equilibrated" without a branch for unequilibrated tissues. So 37 tissues that reached TOTALTIME lost their final state. The user's `~/MainAGDebug.cpp` already writes such states to `_NE.txt`, and their `*_UNFINISHED` parameter files restart from those. The user: "didn't I write a script that does something about this ... Why didn't you use that?"
**Lesson**: If a new driver changes when output is written, every exit path must still save the state, and the user has probably already solved this somewhere. Their Great Lakes home folder (`~/Main*.cpp`, `~/*.sbat`, `~/parameters*`) is the reference for workflow conventions.
**How to Apply**: Before writing or editing a driver, grep `~/Main*.cpp` and the user's batch/parameter files for how each exit path is handled, and copy their convention. Never let a long run exit without writing its final state.


### L007 — Never rename or move a file that a pending or running SLURM job reads
**Mistake**: (2026-10-21) Right after `sbatch --array=1-19 --export=ALL,PARAMS=EQparameters_chunk_resub18.txt`, I renamed the parameter file to fix its name. `batch_EQ_chunked.sbat` reads `$PARAMS` only when each array task starts, so every pending task would have failed. I renamed it back before any task started.
**Lesson**: A batch script's inputs, such as parameter files, binaries and scripts, are read at task start, not at submit time. Choose the final name before `sbatch`.
**How to Apply**: Check the counts and the file name before submitting (e.g. `wc -l` against the intended number). After submitting, leave every file the job touches alone until the job is done.

### L008 — Library sources are CRLF (and some are not pure ASCII): preserve line endings and search with grep -a
**Mistake**: (2026-10-22) I edited `AGPerimTissue.cpp` with Python in text mode. That silently converted the whole file from CRLF to LF, and `git diff` showed 7305 changed lines. The Edit tool also added LF lines to a header. Separately, a plain `grep` skipped `AGCell.cpp` (CRLF plus non-ASCII bytes) as binary, so I briefly concluded that the force code had no QUARTIC2.0 branch.
**Lesson**: Check `file <name>` and the HEAD version's line endings before editing. Edit at the byte level (`open(p,'rb')`, with `\r\n` in the search strings), and afterwards check that `git diff --stat` shows only the intended lines.
**How to Apply**: After any library edit, run `file` and `git diff --stat`. Search library sources with `grep -a`. Don't conclude that something is absent from a search that may have skipped files.

### L009 — Undo a trial move by restoring the saved state (Evolve's backtrack), not by an opposite move; factor shared multipliers out of branches
**Mistake**: (2026-10-24, user code review) `PushAlongMostNegativeModed` undid its trial pushes by pushing back (+s, -2s, +s). That restores the positions only up to round-off, and it left `currentFreeEnergy` at the pre-push value. The user asked: "don't you have to reset those afterwards? UpdateOldMechanicalVariables or something." In the same review, I had put the per-cell `a0 *` into every branch of the pressure-type if-chain instead of applying it once after the chain.
**Lesson**: The codebase already has a save/restore pattern: Evolve saves rold, LXold/LYold and the edge flags, and its backtrack calls `Update_Old_Mechanical_Variables`. Reuse it for any trial move. Keep the energy bookkeeping members consistent with the state that is kept. A factor shared by every branch goes after the branches.
**How to Apply**: Before writing any code that moves vertices, read how Evolve saves and restores state, and copy it. Test that a restore is bit-exact.

### L010 — Check the units of every formula I propose; label options by formulas with consistent units, not bare scalings
**Mistake**: (2026-10-24) I presented push-size options as "1/|λ|" and "1/√|λ|" scalings. The user pointed out that the units don't make sense: s = MINLENGTH/|λ| is not a length. Only s = √(2ΔE/|λ|) is dimensionally consistent. The user also had to ask what E, t_x and t_y were, because I used symbols without defining them.
**Lesson**: Every formula I propose must be dimensionally consistent, written out in full, and come with a definition of every symbol and its units.
**How to Apply**: Before proposing a formula, write its units next to it. Define each variable the first time it appears, in the message itself.

### L011 — When asked to compare shear moduli, give the moduli: a table of G, not a figure or a commentary
**Mistake**: (2026-10-26) Asked "how do MTQ3 and QuTs compare to T2ss?", I made and viewed a figure and led with ratios and mechanism. The user wanted the shear moduli themselves compared.
**Lesson**: A comparison question means numbers first: G mean ± SEM (n seeds) per p0 for each set, side by side. Plots, ratios and explanations come only if asked.
**How to Apply**: Answer with the G table straight from the summary TSV. Don't make or view figures unless the user asks for one.

### L012 — File-transfer commands for the user: explicit full paths, no brace expansion, no extra steps
**Mistake**: (2026-10-26) I gave scp commands with {a,b} brace globs and an unnecessary mkdir. The braces did nothing on the user's setup, and the command was too long.
**Lesson**: The user copies figures to WSL at ~/GreatLakesImages/ (it already exists). Give plain scp commands with every remote file written out in full.
**How to Apply**: One scp per folder, each file as `abmondal@greatlakes.arc-ts.umich.edu:<full path>`, destination `~/GreatLakesImages/`. No braces, no mkdir, no tar.
