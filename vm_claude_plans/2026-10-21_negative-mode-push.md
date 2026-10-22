**Status:** PLANNING STUB, NOT APPROVED (2026-10-21). The user wants this at the top of the todo list and will discuss it before any design is final. Nothing implemented.

# Plan: detect negative Hessian modes, push along the most negative one, re-equilibrate, and record stuck tissues

## Why
A negative Hessian eigenvalue means the "equilibrated" tissue is at a saddle, not a minimum: one direction lowers the energy.
The shear modulus is defined at a minimum. With lambda < 0 the nonaffine term -xi^2/lambda flips sign and RAISES G, so the
value is not a shear modulus. The user: a negative mode is meaningful and must not be ignored.
- Today: the mode loop in `calculateShearModulusWithEigenDecomposition(d)` (AGPerimTissue.cpp ~2926 / ~3180) sends every
  lambda <= threshold, negative ones included, to the floppy test, which returns 0 or drops the mode. Negative modes are never logged.
- Most likely in QZ20 (QUARTIC2.0, C4 = 1/2, A0 = 2 Abar): the per-cell area energy is concave for |a - 1| < 0.816 and the
  mean area is a = 0.5, so cells can lower their energy by trading area.

## User decisions (2026-10-21)
- Zero modes are classified with an absolute bound |lambda| <= threshold. Negative mode: lambda < -threshold.
  The threshold lives in the library, never in a Main file.
- If a negative mode is found: take the MOST negative mode, push the tissue along it, and re-equilibrate.
- The push borrows from the T1 code (see below).
- **If the tissue stays stuck in negative curvature, it must be RECORDED** (the user wants to know).

## The T1 analogue (to discuss)
`Tissue::HelperResolveVertices_ResolveFourfoldVertex` (TissueClass.cpp) resolves a fourfold vertex after a T1 like this:
it tests the candidate resolutions for stability (net force vs. the new edge's tension), picks the stable one, and places the two
new vertices NEWEDGELENGTH/2 on either side of the old vertex, along the force direction (`HelperResolveVertices_UpdatePtrs`,
TissueClass.cpp:2006-2007). Minimization then continues from there.
Proposed analogue for a negative mode u (unit eigenvector, most negative lambda):
1. Amplitude: scale u so the largest single-vertex displacement is a small fixed length (NEWEDGELENGTH-like; value to discuss).
2. Direction: try both +u and -u, and keep the one with the lower energy (the "stable resolution" choice, like T1's stability test).
3. Re-equilibrate with the chunked driver (same CAP/chunk rules), then rerun the Hessian.
4. Repeat up to N_push times (N to discuss).

## Recording (required)
Per tissue, a ledger, e.g. `<name>_negmode.tsv`, with one row per Hessian pass: pass number, number of negative modes, most negative
lambda, lambda_max, |xi . u|, energy before and after the push, chosen sign, push amplitude, and the outcome
(`MIN` = no negative mode, `PUSHED`, `STUCK` = still negative after N_push pushes or the push does not lower the energy).
A per-set summary table reports the fraction of tissues that are MIN, recovered after k pushes, or STUCK, by p0.
STUCK tissues are flagged and kept out of G averages, but their G and lowest eigenvalues are still written out.

## Open questions for the discussion
- Push amplitude and N_push; whether the energy test in step 2 needs a minimum required drop.
- Do the push in the library (a new AGPerimTissue function) or in the driver? Thresholds stay in the library either way.
- Single vs. double precision eigensolver for this check (use the double path).
- Whether T1s during re-equilibration are allowed (they are by default: T1ON = 1).
- Relation to the QUARTIC2.0 plan: QZ20 is the main customer. This plan should be done before the QZ20 Hessian runs.
