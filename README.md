# CeeSectionBuckling

Elastic local and distortional buckling of lipped C (Cee) sections, gross and net (web hole) models, by finite strip analysis of the actual rounded-corner centerline geometry, with the critical points located by the buckling mode identification procedure of the RunToSolve ESR-1538 report (Scottsdale Construction Systems, 2026), implemented with [BucklingModeIdentification.jl](https://github.com/runtosolve/BucklingModeIdentification.jl).

```julia
using CeeSectionBuckling
material   = CeeSectionBuckling.Material(29500.0, 0.30)                 # E (ksi), ν
dimensions = CeeSectionBuckling.Dimensions(0.0713, 1.0, 4.0, 3.62, 0.1070)   # t, lip L, flange B, web H (out-to-out), inside radius r
Mcrl = CeeSectionBuckling.calculate_Mcrℓ_xx(dimensions, material)           # Section; Mcrl.results.Rcr (kip-in.), .Lcr (in.)
println(CeeSectionBuckling.describe(Mcrl.results))
props = CeeSectionBuckling.calculate_properties(dimensions, material)       # section properties + Pcrℓ, Pcrd, Mcrd_xx, Mcrℓ_xx, Mcrd_yy_pos, Mcrℓ_yy_pos, Mcrℓ_yy_neg
Mcrl_net = CeeSectionBuckling.calculate_Mcrℓ_xx_net_section(dimensions, material, 1.5, 4.0)   # 1.5 in. wide x 4.0 in. long web punchout, zero-thickness net section
Mcrd_net = CeeSectionBuckling.calculate_Mcrd_xx_net_section(dimensions, material, 4.0)   # 4.0 in. long web punchout, reduced web thickness (AISI S100 App. 2)
```

## How the critical points are located

Every value is the conventional (unconstrained) finite strip signature curve of the rounded-corner geometry; the constrained finite strip method (cFSM) only ever decides *where* on that curve to read.

**Local buckling** (`identify_local_buckling`)
1. Conventional sweep over [0.5W, 2W], W = widest flat element of the section (web or flange, between corner arcs), 31 log-spaced half-wavelengths with the range ends sampled exactly (plus one guard point outside each end).
2. The first true trough inside the range, R[i-1] > R[i] < R[i+1], whose buckled shape is L-dominated (largest G/D/L/O participation is L) is the local mode (`method = :signature_trough`). Values at the range ends do not matter. For deep webs with narrow flanges the distortional trough itself falls inside [0.5W, 2W]; the participation keeps it from being read as the local mode.
3. If the curve has no L-dominated trough in the range (local and distortional modes interact and the curve only falls through the range), the pure-local cFSM curve of the same rounded geometry (elastic-corner model) gives Lcrℓ at its minimum (the constrained sweep is extended beyond the window when its minimum sits on an end) and the reported value is the signature-curve value at that Lcrℓ (`method = :cFSM_pure_mode`).

**Distortional buckling** (`identify_distortional_buckling`)
1. Closed-form AISI S100 Appendix 2 half-wavelength estimate (`calculate_Lcrd`) seeds the search.
2. Conventional sweeps: 5 half-wavelengths over [0.5, 2] × estimate, then 10 over [0.7, 1.3] × that minimum. The minimum of the merged curve is Lcrd when it is a true trough and D-dominated (`:signature_trough`). For deep webs the window contains the web's local trough, which is lower than the distortional one and is rejected by the participation check.
3. Otherwise a wide sweep over [0.2, 5] × estimate is decomposed into G/D/L/O participations and the D-dominated trough is refined; failing that (no distortional trough at all, e.g. in compression where global buckling undercuts it), the pure-distortional cFSM minimum is used (`:cFSM_pure_mode`).

**Web holes.** Local (`identify_local_buckling_at_hole`): one zero-thickness web element across the hole width (Moen & Schafer); Lcrℓ = min(Lh, Lcrℓ,cFSM), where Lh is the hole length and Lcrℓ,cFSM is the pure-local cFSM minimum of the net-section geometry, and the reported value is the conventional signature-curve value at that Lcrℓ (`method = :hole_length` or `:cFSM_pure_mode`). The conventional sweep over the gross section's [0.5W, 2W] is recorded for figures only. Distortional: reduced web thickness tr = t(1 − Lh/Lcrd)^(1/3) over the whole web, evaluated at the gross section's identified Lcrd (`method = :specified_length`).

`results.identification` (an `Identification`) records the sweep, its troughs, the cFSM pure-mode curve and minimum, the conventional value at the cFSM half-wavelength, and the G/D/L/O participation at the reported point, so a figure like the reference report's Figure 2 can be drawn from it. `results.model` is the CUFSM model evaluated at the reported half-wavelength, so `results.model.shapes[1]` is the reported buckled shape.

## Tests

`test/mode_identification_test.jl` checks the procedure on the MRI 362S400-68 stud (gross and net-section modes); the other `test/*.jl` scripts are usage examples with mode-shape figures.
