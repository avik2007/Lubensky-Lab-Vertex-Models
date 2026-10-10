# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Workflow Orchestration

### 1. Plan Mode Default
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions).
- If something goes sideways, STOP and re-plan immediately.
- Use plan mode for verification steps, not just building.
- Write detailed specs upfront to reduce ambiguity.

### 2. Subagent Strategy
- Use subagents liberally to keep main context window clean.
- Offload research, exploration, and parallel analysis to subagents.
- For complex problems, throw more compute at it via subagents.
- One task per subagent for focused execution.

### 3. Self-Improvement Loop
- After ANY correction from the user: update `claude_actions/claude_lessons.md` with the pattern.
- Write rules for yourself that prevent the same mistake.
- Ruthlessly iterate on these lessons until mistake rate drops.
- Review `claude_actions/claude_lessons.md` at session start for relevant context.
- When lessons accumulate, update this CLAUDE.md accordingly.

### 4. Verification Before Done
- Never mark a task complete without proving it works.
- Diff behavior between main and your changes when relevant.
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness.

### 5. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution."
- Skip this for simple, obvious fixes — don't over-engineer.
- Challenge your own work before presenting it.

### 6. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding.
- Point at logs, errors, failing tests — then resolve them.
- Zero context switching required from the user.
- Go fix failing CI tests without being told how.

---

## Task Management

- **Plan First**: Write plan to `claude_actions/claude_todo.md` with checkable items.
- **Verify Plan**: Check in before starting implementation.
- **Track Progress**: Mark items complete as you go.
- **Explain Changes**: High-level summary at each step.
- **Document Results**: Add review section to `claude_actions/claude_todo.md`.
- **Capture Lessons**: Update `claude_actions/claude_lessons.md` after corrections.

---

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **Verbose Comments**: Every new piece of code must have verbose comments explaining what it does, why it exists, and how it works — clear enough for a human reader unfamiliar with the codebase to understand.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Only touch what's necessary. No side effects with new bugs.

---

## Project Context

### Scientific Workflow
Scientific planning and brainstorming happen in a separate LLM. Prompts are generated there and fed here for implementation. Claude Code's role is **execution, not scientific ideation** — translate well-specified prompts into working code and analysis.

### Branch Scope: `MondalShearDev`
The primary purpose of this branch is to **fix and improve shear modulus and Hessian calculations**. Relevant code lives in:
- `AG Vertex Model/AGVertexModel/AGStiffness/` — shear/elastic stiffness
- `AG Vertex Model/AGVertexModel/AGPerimModuli/` — elastic moduli
- `AG Vertex Model/AGVertexModel/AGPerimHessian_test/` — Hessian extraction
- `Python Code Folder/script_TO_CALC_SM_FROM_CPP.py` — shear modulus post-processing

**Do not touch code outside these areas unless explicitly asked.** Changes to unrelated simulation variants, MATLAB generators, or analysis utilities are out of scope.

---

## Project Overview

2D vertex model simulations of epithelial tissue dynamics (Lubensky group, University of Michigan). MATLAB generates initial conditions/input files; C++ runs simulations; Python analyzes results.

**Full file structure, class hierarchy, data flow, and environment setup:** see [`vm_file_structure.md`](vm_file_structure.md).

### File Structure Maintenance
**Whenever you add or remove files or directories**, update `vm_file_structure.md` immediately to reflect the change. This is the single source of truth for repo layout.
