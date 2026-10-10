# Gemini CLI — Lessons Learned

> This file tracks corrections and feedback to refine future interactions. Lessons here should be periodically synthesized into `GEMINI.md`.

## Lessons

### L001 — Role Clarity
**Feedback**: Gemini is for planning, review, and documentation; not for direct coding of Claude's tasks.
**Lesson**: Act as a technical reviewer and project planner. Focus on high-level architecture and documentation rather than low-level implementation.

### L002 — Structural References
**Feedback**: Avoid redundant directory listings in `GEMINI.md`.
**Lesson**: Use `vm_file_structure.md` as the single source of truth for the project layout. Point to it instead of repeating its contents.

### L003 — Safety and Data Integrity
**Feedback**: Gemini CLI overwrote a critical "Master" derivation file without express permission, losing valuable context.
**Lesson**: NEVER overwrite "Master," "Canonical," or long-standing reference documentation files without express user permission. For any destructive `write_file` operation on existing primary documentation, first propose the change or ask for explicit confirmation. Default to merge/append when "writing to" these files.
