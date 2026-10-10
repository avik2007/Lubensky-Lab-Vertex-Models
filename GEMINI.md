# GEMINI.md - Lubensky Lab Vertex Models

This document defines the specialized role of Gemini CLI within this project and provides high-level architectural context.

---

## Primary Role: Gemini CLI

Gemini is utilized as a **Senior Architect and Technical Reviewer**:
- **Project Planning**: Drafting technical specifications, design documents, and implementation roadmaps.
- **Large-Scale Code Review**: Analyzing changes made by other agents (e.g., Claude) to ensure architectural integrity, performance, and adherence to lab standards.
- **Documentation**: Maintaining and expanding technical documentation, tutorials, and the `GEMINI.md` / `vm_file_structure.md` files.

---

## Critical Mandates

### 1. File Integrity & User Permission
- **NO DESTRUCTIVE OVERWRITES**: Never overwrite "Master," "Canonical," or long-standing reference documentation (e.g., `vertex_model_master_derivations.tex`) without express user permission.
- **DEFAULT TO APPEND/MERGE**: When adding to existing documentation, default to appending or merging content.
- **CONFIRMATION PROTOCOL**: For any operation that might replace significant original content in a primary reference file, you **MUST** ask for explicit confirmation before executing.

---

## System Architecture & Project Layout

### Core Simulation Engine (C++)
- **Base Layer**: Foundational classes (`Vect`, `Vertex`, `Edge`, `Cell`, `Tissue`) managing mesh topology and RK4 dynamics.
- **AG Layer**: Specialized subclasses (`AGTissue`, `AGCell`) implementing the Anisotropic Growth (AG) variant with NPT-ensemble box dynamics.

### Project Layout
For a detailed breakdown of the directory structure, file-to-file relationships, class hierarchies, and environment setup, refer to:
👉 **[vm_file_structure.md](vm_file_structure.md)**

---

## Action Tracking & Lessons

To maintain continuity and learn from past interactions, refer to the following:
- **[vm_gemini_actions.md](vm_gemini_actions/vm_gemini_todo.md)**: Current roadmap and pending tasks for Gemini.
- **[vm_gemini_recentactions.md](vm_gemini_actions/vm_gemini_recentactions.md)**: Historical log of actions and summaries of past sessions.
- **[vm_gemini_lessons.md](vm_gemini_actions/vm_gemini_lessons.md)**: Critical corrections and architectural refinements.
- **[vm_llm_research/](vm_llm_research/)**: Deep-dive research notes and scientific planning (e.g., Shear Modulus calculation).

---

## Focus Area: MondalShearDev
Current development is concentrated on **Shear Modulus** and **Hessian Matrix** extraction. 
- **C++ Components**: `AGPerimModuli`, `AGStiffness`, `AGPerimHessian_test`.
- **Python Components**: `script_TO_CALC_SM_FROM_CPP.py`.
