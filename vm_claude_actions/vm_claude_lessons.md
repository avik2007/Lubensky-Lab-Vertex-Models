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
