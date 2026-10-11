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
