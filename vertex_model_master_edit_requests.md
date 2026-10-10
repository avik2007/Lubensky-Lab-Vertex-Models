# vertex_model_master_edit_requests.md

This file tracks requested improvements and clarifications for `vertex_model_master_derivations.tex`.

...

## 7. Diagram Legibility and Logic Refinement
**Status:** Completed.

## 8. Pivot to Standard Linear Response (Section A.6)
**Request:** Replace the "Extended Hessian" (Section A.7) implementation guide with the "Standard Linear Response" (Section A.6) method. 
**Constraints:** 
- $\gamma$ must be an external parameter, not a degree of freedom.
- Use the spatial Hessian $\mathcal{H}$ ($DN_v \times DN_v$) exclusively for diagonalization.
- Implement the non-affine term via eigenvector projection: $\sum \frac{1}{\omega^2} (\vb{\Xi} \cdot \vb{u})^2$.
- Include explicit zero-mode filtering and floppy-mode detection.
**Status:** Pending.
