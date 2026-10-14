module CeeSectionBuckling


using CrossSectionGeometry, CUFSM, AISIS100, SectionProperties, BucklingModeIdentification

const BMI = BucklingModeIdentification


struct Material

    E
    ν

end

struct Dimensions

    t

    L   #lip
    B   #flange
    H   #web
    r   #inside radius

end


struct Load

    P
    Mxx
    Mzz
    M11
    M22

end


#How a critical local or distortional buckling point was located on the conventional
#(unconstrained) finite strip signature curve of the actual rounded-corner geometry,
#following the RunToSolve ESR-1538 (Scottsdale, 2026) report procedure:
#
#  local (mode = :local)
#    1. Conventional sweep over [0.5W, 2W], W = widest flat element of the section (web
#       or flange, between corner arcs), 31 log-spaced half-wavelengths with 0.5W and 2W
#       sampled exactly, plus one guard point just outside each end so that a trough at
#       a range end can be recognised.
#    2. Rule: the FIRST true local minimum (trough) inside the range, R[i-1] > R[i] < R[i+1],
#       is the local buckling mode -> method = :signature_trough. The values at the range
#       ends do not matter (after the trough the curve may fall again towards the
#       distortional / global modes).
#    3. Otherwise the curve declines monotonically through the range because the local
#       and distortional modes interact: the pure-local constrained finite strip (cFSM)
#       curve of the SAME rounded geometry (BucklingModeIdentification.jl, elastic-corner
#       model) locates Lcrℓ at its minimum and the reported value is the conventional
#       signature-curve value at that Lcrℓ (always at or below the constrained minimum)
#       -> method = :cFSM_pure_mode.
#
#  distortional (mode = :distortional)
#    1. Closed-form AISI S100 Appendix 2 half-wavelength estimate Lcrd,est (Lcr_estimate).
#    2. Conventional sweep, 5 half-wavelengths over [0.5, 2.0] x Lcrd,est, then 10 over
#       [0.7, 1.3] x that minimum; Lcrd is the minimum of the merged curve when it is a
#       true trough -> method = :signature_trough.
#    3. Otherwise a wide sweep over [0.2, 5] x Lcrd,est is decomposed into G/D/L/O
#       participations and the D-dominated trough with the lowest load factor is refined
#       (method = :signature_trough); failing that, Lcrd is the pure-distortional cFSM
#       minimum and the reported value is the conventional curve there
#       (method = :cFSM_pure_mode).
#
#  local at a web hole (identify_local_buckling_at_hole, net-section models)
#    Lcrℓ = min(Lh, Lcrℓ,cFSM): the pure-local cFSM minimum of the net-section geometry,
#    limited to the hole length Lh because the local buckling half-wave of the net
#    section cannot be longer than the hole. The reported value is the conventional
#    signature-curve value at that Lcrℓ -> method = :hole_length when Lh governs,
#    :cFSM_pure_mode otherwise.
#
#  method = :specified_length: the value was evaluated at a half-wavelength fixed by
#  another analysis (the reduced-thickness hole models are evaluated at the gross
#  section's Lcrd).
#
#The cFSM pure-mode curve (curve_cFSM at lengths_cFSM, minimum Lcr_cFSM / Rcr_cFSM, and
#the conventional value at that half-wavelength, Rcr_at_Lcr_cFSM) is always recorded so
#that a figure like the reference report's Figure 2 can be drawn from this record alone.
#`participation` is the G/D/L/O modal participation (percent) of the conventional
#buckled shape at the reported point.
struct Identification

    mode::Symbol            # :local or :distortional
    method::Symbol          # :signature_trough, :cFSM_pure_mode, :hole_length, or :specified_length
    W                       # widest flat element (local search-range basis), nothing for distortional
    Lcr_estimate            # closed-form Lcrd estimate (distortional), nothing for local
    search_range            # (Lmin, Lmax) of the conventional sweep
    lengths                 # conventional sweep half-wavelengths (guard points included)
    curve                   # conventional load factors at `lengths`
    sweep_model             # CUFSM.Model of the conventional sweep (buckled shapes at every length)
    troughs                 # indices into `lengths` of every pointwise trough R[i-1] > R[i] < R[i+1]
    cFSM_subspace           # :L or :D
    lengths_cFSM            # half-wavelengths of the pure-mode (cFSM) curve, nothing if not run
    curve_cFSM              # pure-mode load factors (NaN where no positive load factor exists)
    Lcr_cFSM                # refined pure-mode minimum half-wavelength, nothing if not run / not found
    Rcr_cFSM                # pure-mode minimum load factor
    Rcr_at_Lcr_cFSM         # conventional signature-curve value at Lcr_cFSM
    participation           # (G, D, L, O) percent of the conventional mode at the reported point

end


#`model` is the CUFSM.Model evaluated AT the reported half-wavelength (model.lengths ==
#[Lcr], get_load_factor(model, 1) == [Rcr]), so model.shapes[1] is the reported buckled
#shape. `identification` documents how (Lcr, Rcr) was located (an `Identification`), or
#is `nothing` for the low-level sweep functions that simply return the minimum of the
#half-wavelengths they were given.
struct Results

    model
    Lcr
    Rcr
    identification

end

Results(model, Lcr, Rcr) = Results(model, Lcr, Rcr, nothing)


struct Section

    label
    material
    dimensions
    load
    results

end




function get_section_coordinates(dimensions)

    (;
    t,

    L,
    B,
    H,
    r  #inside radius
    ) = dimensions

    section_dimensions = [L, B, H, B, L]
    r = [r+t, r+t, r+t, r+t]
    n = [3, 3, 3, 3, 3]
    n_r = [3, 3, 3, 3];
    n_r = [7, 7, 7, 7];
    θ = [π/2, π, -π/2, 0.0, π/2]

    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, r, n_r, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3] - section_dimensions[1]))

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    coordinates = (X=X, Y=Y)

    return coordinates

end



function get_straight_corner_section_coordinates(dimensions)

    (;
    t,

    L,
    B,
    H,
    r  #inside radius
    ) = dimensions

    section_dimensions = [L, B, H, B, L]
    n = [1, 1, 1, 1, 1]
    θ = [π/2, π, -π/2, 0.0, π/2]

    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3] - section_dimensions[1]))

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    coordinates = (X=X, Y=Y, L=section_dimensions)

    return coordinates

end


function calculate_section_properties(coordinates, t)

    X = coordinates.X
    Y = coordinates.Y

    t_all = fill(t, length(X)-1)

    section_properties = SectionProperties.open_thin_walled(X, Y, t_all)

    return section_properties

end



# ---------------------------------------------------------------------------
# Conventional finite strip sweeps (low level): minimum over the given half-wavelengths
# ---------------------------------------------------------------------------

#Conventional (unconstrained) CUFSM analysis of the centerline `X, Y` with element
#thickness `t` (scalar, or one value per element) at the given half-wavelengths.
function conventional_sweep(X, Y, t, loads, material, lengths)

    (;E,
    ν
    ) = material

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    constraints = []
    springs = []

    neigs = 1

    model = CUFSM.Tools.open_section_analysis(X, Y, t, collect(Float64, lengths), E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    return model

end


function calculate_buckling_properties(coordinates, dimensions, loads, material, lengths)

    t = dimensions.t

    model = conventional_sweep(coordinates.X, coordinates.Y, t, loads, material, lengths)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)
    Rcr = minimum(Rcr_curve)

    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results

end


#Same as `calculate_buckling_properties`, but takes the element thicknesses directly
#(scalar or per-element vector) instead of pulling a single thickness from
#`dimensions`, so it can be used with the net section's zeroed-out element
#thicknesses at a web hole or the reduced web thickness of the distortional hole model.
function calculate_net_section_buckling_properties(coordinates, t, loads, material, lengths)

    model = conventional_sweep(coordinates.X, coordinates.Y, t, loads, material, lengths)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)
    Rcr = minimum(Rcr_curve)

    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results

end



# ---------------------------------------------------------------------------
# Buckling mode identification on the actual rounded-corner geometry
# ---------------------------------------------------------------------------

#Widest flat element of the centerline (web or flange, measured between corner arcs),
#the W that sets the local buckling search range [0.5W, 2W].
function max_flat_width(X, Y)

    runs = BMI.straight_runs(X, Y)

    return maximum(sum(hypot.(diff(X[r.start:r.stop+1]), diff(Y[r.start:r.stop+1]))) for r in runs)

end


#Indices of every pointwise trough of a curve, R[i-1] > R[i] < R[i+1].
troughs(R) = [i for i in 2:(length(R) - 1) if R[i-1] > R[i] < R[i+1]]


#CUFSM arrays (prop, node, elem) of the same geometry, thickness, and reference stress
#for BucklingModeIdentification's constrained (cFSM) analyses and decompositions.
function cfsm_model(X, Y, t, loads, material)

    return BMI.cufsm_model(X, Y, t, material.E; ν = material.ν, P = loads.P, Mxx = loads.Mxx, Mzz = loads.Mzz, M11 = loads.M11, M22 = loads.M22)

end


#Pure-mode (constrained finite strip) curve of `subspace` (:L or :D) at `lengths`, and
#its minimum refined with `n_refine` log-spaced half-wavelengths between the neighbours
#of the coarse minimum. While the minimum sits on an end point of the curve the sweep is
#extended a decade beyond that end (at most `max_extensions` times), because the
#constrained mode's minimum can lie outside the conventional search window - e.g. the
#flange-tip local mode of a very deep, narrow-flange section in weak-axis bending has a
#half-wavelength set by the flange, well below 0.5W of the web. Lcr/Rcr are `nothing`
#when the constrained problem has no positive load factor anywhere.
function pure_mode_minimum(node, elem, prop, lengths, subspace; n_refine = 15, n_extend = 12, max_extensions = 3)

    lengths = collect(Float64, lengths)
    curve = BMI.pure_mode_curve(node, elem, prop, lengths, [subspace])

    finite = findall(isfinite, curve)
    isempty(finite) && return (lengths = lengths, curve = curve, Lcr = nothing, Rcr = nothing)

    i = finite[argmin(curve[finite])]

    extensions = 0
    while (i == 1 || i == length(lengths)) && extensions < max_extensions
        L_extend = i == 1 ? exp.(range(log(lengths[1] / 10), log(lengths[1]), length = n_extend + 1))[1:end-1] :
                            exp.(range(log(lengths[end]), log(10 * lengths[end]), length = n_extend + 1))[2:end]
        curve_extend = BMI.pure_mode_curve(node, elem, prop, L_extend, [subspace])
        order = sortperm(vcat(lengths, L_extend))
        lengths = vcat(lengths, L_extend)[order]
        curve = vcat(curve, curve_extend)[order]
        finite = findall(isfinite, curve)
        i = finite[argmin(curve[finite])]
        extensions += 1
    end

    L_refine = exp.(range(log(lengths[max(i - 1, 1)]), log(lengths[min(i + 1, end)]), length = n_refine))
    curve_refine = BMI.pure_mode_curve(node, elem, prop, L_refine, [subspace])

    finite_refine = findall(isfinite, curve_refine)
    isempty(finite_refine) && return (lengths = lengths, curve = curve, Lcr = lengths[i], Rcr = curve[i])

    j = finite_refine[argmin(curve_refine[finite_refine])]

    return (lengths = lengths, curve = curve, Lcr = L_refine[j], Rcr = curve_refine[j])

end


#G/D/L/O modal participation (percent) of the conventional buckled shape at
#half-wavelength `L`.
function participation_at(node, elem, prop, L)

    decomposition = BMI.decompose_modes(node, elem, prop, [L])
    p = decomposition.participation[1, :]

    return (G = p[1], D = p[2], L = p[3], O = p[4])

end


#Class with the largest participation, :G, :D, :L, or :O.
dominant(p) = (:G, :D, :L, :O)[argmax((p.G, p.D, p.L, p.O))]


#Elastic local buckling of the centerline `X, Y` with element thickness `t` (scalar or
#per element) under `loads`: first L-dominated trough of the conventional signature
#curve inside [0.5W, 2W], otherwise the conventional value at the cFSM pure-local
#half-wavelength. See `Identification` for the procedure. `W` defaults to the widest
#flat of `X, Y`; the net-section models pass the gross section's W so that gross and
#net models are searched over the same range. `cfsm = false` skips the constrained
#analysis (then a curve with no L-dominated trough in the range is an error).
function identify_local_buckling(X, Y, t, loads, material; W = nothing, range_factors = (0.5, 2.0), n_range = 31, guard = 1.05, cfsm = true)

    W === nothing && (W = max_flat_width(X, Y))
    Lmin, Lmax = range_factors .* W

    inner = exp.(range(log(Lmin), log(Lmax), length = n_range))
    lengths = vcat(Lmin / guard, inner, Lmax * guard)   #guard points: a trough at 0.5W or 2W needs a neighbour beyond it

    sweep = conventional_sweep(X, Y, t, loads, material, lengths)
    R = CUFSM.Tools.get_load_factor(sweep, 1)
    trough_indices = troughs(R)   #never the guard points (a trough needs neighbours on both sides)

    prop, node, elem = cfsm_model(X, Y, t, loads, material)

    #a trough is the local mode only if the buckled shape there is L-dominated: for deep
    #webs with narrow flanges the distortional trough falls inside [0.5W, 2W] and would
    #otherwise be mistaken for the local mode (Mcrℓ = Mcrd)
    local_troughs = Int[]
    participation = nothing
    for i in trough_indices
        p = participation_at(node, elem, prop, lengths[i])
        if dominant(p) == :L
            push!(local_troughs, i)
            participation === nothing && (participation = p)
        end
    end

    if cfsm
        pm = pure_mode_minimum(node, elem, prop, lengths, :L)
        Rcr_at_Lcr_cFSM = pm.Lcr === nothing ? nothing : CUFSM.Tools.get_load_factor(conventional_sweep(X, Y, t, loads, material, [pm.Lcr]), 1)[1]
    else
        pm = (lengths = nothing, curve = nothing, Lcr = nothing, Rcr = nothing)
        Rcr_at_Lcr_cFSM = nothing
    end

    if !isempty(local_troughs)
        i = first(local_troughs)
        Lcr = lengths[i]
        method = :signature_trough
    elseif pm.Lcr !== nothing
        Lcr = pm.Lcr
        method = :cFSM_pure_mode
        participation = participation_at(node, elem, prop, Lcr)
    else
        error("local buckling could not be identified: the signature curve has no L-dominated trough in [$(Lmin), $(Lmax)] and no cFSM pure-local minimum is available")
    end

    model = conventional_sweep(X, Y, t, loads, material, [Lcr])
    Rcr = CUFSM.Tools.get_load_factor(model, 1)[1]

    identification = Identification(:local, method, W, nothing, (Lmin, Lmax), lengths, R, sweep, local_troughs, :L,
                                    pm.lengths, pm.curve, pm.Lcr, pm.Rcr, Rcr_at_Lcr_cFSM, participation)

    return Results(model, Lcr, Rcr, identification)

end


#Elastic distortional buckling of the centerline `X, Y` (thickness `t`, scalar or per
#element) under `loads`, searched around the closed-form half-wavelength estimate
#`Lcrd_estimate`: coarse 5-point then fine 10-point conventional sweeps; the minimum must
#be a true trough, otherwise the D-dominated trough of a wide decomposed sweep, otherwise
#the cFSM pure-distortional minimum. See `Identification` for the procedure.
function identify_distortional_buckling(X, Y, t, loads, material, Lcrd_estimate;
                                        coarse = (0.5, 2.0, 5), fine = (0.7, 1.3, 10), wide = (0.2, 5.0, 25), cfsm = true)

    #stage 1: coarse sweep about the estimate
    L_coarse = collect(range(coarse[1] * Lcrd_estimate, coarse[2] * Lcrd_estimate, coarse[3]))
    R_coarse = CUFSM.Tools.get_load_factor(conventional_sweep(X, Y, t, loads, material, L_coarse), 1)
    L1 = L_coarse[argmin(R_coarse)]

    #stage 2: fine sweep about the coarse minimum, merged with the coarse points
    L_fine = collect(range(fine[1] * L1, fine[2] * L1, fine[3]))
    lengths = sort(unique(vcat(L_coarse, L_fine)))
    sweep = conventional_sweep(X, Y, t, loads, material, lengths)
    R = CUFSM.Tools.get_load_factor(sweep, 1)
    trough_indices = troughs(R)

    i_min = argmin(R)
    method = :signature_trough
    Lcr = lengths[i_min]

    prop, node, elem = cfsm_model(X, Y, t, loads, material)

    #the minimum must be a true trough AND the buckled shape there must be D-dominated:
    #for deep webs with narrow flanges the window about the Appendix 2 estimate contains
    #the web's local trough (at about W, below the estimate), which is lower than the
    #distortional one and would otherwise be taken as Lcrd
    participation = i_min in trough_indices ? participation_at(node, elem, prop, Lcr) : nothing
    accepted = participation !== nothing && dominant(participation) == :D

    if !accepted

        #stage 3: wide sweep decomposed into G/D/L/O at every half-wavelength, so the
        #distortional trough can be told from the local and global ones
        L_wide = exp.(range(log(wide[1] * Lcrd_estimate), log(wide[2] * Lcrd_estimate), length = wide[3]))
        decomposition = BMI.decompose_modes(node, elem, prop, L_wide)
        minima = BMI.characteristic_minima(decomposition)

        if minima.distortional_min !== nothing
            Ld = L_wide[minima.distortional_min]
            L_fine = collect(range(fine[1] * Ld, fine[2] * Ld, fine[3]))
            lengths = sort(unique(vcat(L_wide, L_fine)))
            sweep = conventional_sweep(X, Y, t, loads, material, lengths)
            R = CUFSM.Tools.get_load_factor(sweep, 1)
            trough_indices = troughs(R)
            candidates = filter(i -> fine[1] * Ld * (1 - 1e-9) <= lengths[i] <= fine[2] * Ld * (1 + 1e-9) &&
                                     dominant(participation_at(node, elem, prop, lengths[i])) == :D, trough_indices)
            #the refined D-dominated trough nearest the decomposed one; the decomposed trough
            #itself (a point of the merged sweep) if the finer grid shows none in that window
            Lcr = isempty(candidates) ? Ld : lengths[candidates[argmin(R[candidates])]]
        else
            method = :cFSM_pure_mode   #resolved below from the pure-distortional curve
        end

    end

    if cfsm || method == :cFSM_pure_mode
        pm = pure_mode_minimum(node, elem, prop, lengths, :D)
        Rcr_at_Lcr_cFSM = pm.Lcr === nothing ? nothing : CUFSM.Tools.get_load_factor(conventional_sweep(X, Y, t, loads, material, [pm.Lcr]), 1)[1]
    else
        pm = (lengths = nothing, curve = nothing, Lcr = nothing, Rcr = nothing)
        Rcr_at_Lcr_cFSM = nothing
    end

    if method == :cFSM_pure_mode
        pm.Lcr === nothing && error("distortional buckling could not be identified: no D-dominated trough of the signature curve about Lcrd = $(Lcrd_estimate) and no cFSM pure-distortional minimum is available")
        Lcr = pm.Lcr
    end

    model = conventional_sweep(X, Y, t, loads, material, [Lcr])
    Rcr = CUFSM.Tools.get_load_factor(model, 1)[1]

    participation = participation_at(node, elem, prop, Lcr)

    identification = Identification(:distortional, method, nothing, Lcrd_estimate, (minimum(lengths), maximum(lengths)), lengths, R, sweep, trough_indices, :D,
                                    pm.lengths, pm.curve, pm.Lcr, pm.Rcr, Rcr_at_Lcr_cFSM, participation)

    return Results(model, Lcr, Rcr, identification)

end


#Elastic local buckling of a net section through a web hole of length `Lh` (along the
#member axis), centerline `X, Y` with element thickness `t` (the zero-thickness hole
#elements of calculate_net_section_element_thicknesses): Lcrℓ = min(Lh, Lcrℓ,cFSM), with
#Lcrℓ,cFSM the minimum of the pure-local constrained finite strip curve of the same
#net-section geometry, and the reported value is the conventional (unconstrained)
#signature-curve value at that Lcrℓ. The conventional sweep over [0.5W, 2W] is recorded
#for figures but does not locate the point. See `Identification`.
function identify_local_buckling_at_hole(X, Y, t, loads, material, Lh; W = nothing, range_factors = (0.5, 2.0), n_range = 31, guard = 1.05)

    W === nothing && (W = max_flat_width(X, Y))
    Lmin, Lmax = range_factors .* W

    inner = exp.(range(log(Lmin), log(Lmax), length = n_range))
    lengths = vcat(Lmin / guard, inner, Lmax * guard)

    sweep = conventional_sweep(X, Y, t, loads, material, lengths)
    R = CUFSM.Tools.get_load_factor(sweep, 1)

    prop, node, elem = cfsm_model(X, Y, t, loads, material)
    pm = pure_mode_minimum(node, elem, prop, lengths, :L)
    pm.Lcr === nothing && error("net-section local buckling could not be identified: no cFSM pure-local minimum is available")
    Rcr_at_Lcr_cFSM = CUFSM.Tools.get_load_factor(conventional_sweep(X, Y, t, loads, material, [pm.Lcr]), 1)[1]

    if Lh < pm.Lcr
        Lcr = Lh
        method = :hole_length
    else
        Lcr = pm.Lcr
        method = :cFSM_pure_mode
    end

    model = conventional_sweep(X, Y, t, loads, material, [Lcr])
    Rcr = CUFSM.Tools.get_load_factor(model, 1)[1]

    participation = participation_at(node, elem, prop, Lcr)

    identification = Identification(:local, method, W, nothing, (Lmin, Lmax), lengths, R, sweep, troughs(R), :L,
                                    pm.lengths, pm.curve, pm.Lcr, pm.Rcr, Rcr_at_Lcr_cFSM, participation)

    return Results(model, Lcr, Rcr, identification)

end


#Conventional value of the centerline `X, Y` (thickness `t`, scalar or per element) at a
#half-wavelength `Lcr` fixed by another analysis (method = :specified_length), with the
#G/D/L/O participation of the mode there.
function evaluate_at_length(X, Y, t, loads, material, Lcr; mode = :distortional, Lcr_estimate = nothing)

    model = conventional_sweep(X, Y, t, loads, material, [Lcr])
    Rcr = CUFSM.Tools.get_load_factor(model, 1)[1]

    prop, node, elem = cfsm_model(X, Y, t, loads, material)
    participation = participation_at(node, elem, prop, Lcr)

    identification = Identification(mode, :specified_length, nothing, Lcr_estimate, (Lcr, Lcr), [Lcr], [Rcr], model, Int[], mode == :local ? :L : :D,
                                    nothing, nothing, nothing, nothing, nothing, participation)

    return Results(model, Lcr, Rcr, identification)

end



# ---------------------------------------------------------------------------
# Closed-form distortional buckling half-wavelength estimate (search seed)
# ---------------------------------------------------------------------------

function calculate_Lcrd(dimensions, material, load_type)

    (;t,
    L,
    B,
    H,
    r ) = dimensions

    (;E,
    ν) = material

    CorZ = 0
    θ_top = 90.0
    #Calculate top flange + lip section properties.
    CorZ = 0

    b = B - t
    d = L - t/2
    θ = 90.0
    ho = H
    μ = ν
    E = 29500.0
    G = E / (2 * (1 + μ))
    kϕ = 0.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)

    #Calculate the purlin distortional buckling half-wavelength.
    Lm = 999999999.0

    if load_type == "P"
        Lcrd = AISIS100.v16S3.appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    elseif load_type == "M"
        Lcrd = AISIS100.v16S3.appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    end

    return Lcrd

end



# ---------------------------------------------------------------------------
# Gross-section modes
# ---------------------------------------------------------------------------

#Local ("L") or distortional ("D") elastic buckling of the gross section under `load`.
#Distortional searches are seeded with the closed-form Lcrd estimate for compression
#(load.P != 0) or bending (otherwise).
function calculate_Rcr(label, load, mode_type, dimensions, material)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates

    if mode_type == "L"
        results = identify_local_buckling(coordinates.X, coordinates.Y, dimensions.t, load, material)
    elseif mode_type == "D"
        load_type = load.P != 0 ? "P" : "M"
        Lcrd_estimate = calculate_Lcrd(dimensions, material, load_type)
        results = identify_distortional_buckling(coordinates.X, coordinates.Y, dimensions.t, load, material, Lcrd_estimate)
    else
        error("mode_type must be \"L\" or \"D\", got $(mode_type)")
    end

    section = Section(label, material, dimensions, load, results)

    return section

end


function calculate_Pcrℓ(dimensions, material)

    return calculate_Rcr("Pcrℓ", Load(1.0, 0.0, 0.0, 0.0, 0.0), "L", dimensions, material)

end


function calculate_Pcrd(dimensions, material)

    return calculate_Rcr("Pcrd", Load(1.0, 0.0, 0.0, 0.0, 0.0), "D", dimensions, material)

end


function calculate_Mcrd_xx(dimensions, material)

    return calculate_Rcr("Mcrd_xx", Load(0.0, 1.0, 0.0, 0.0, 0.0), "D", dimensions, material)

end


function calculate_Mcrℓ_xx(dimensions, material)

    return calculate_Rcr("Mcrℓ_xx", Load(0.0, 1.0, 0.0, 0.0, 0.0), "L", dimensions, material)

end


#this is not the traditional strong axis distortional buckling, still important though in some cases for weak axis bending
function calculate_Mcrd_yy_pos(dimensions, material)

    return calculate_Rcr("Mcrd_yy_pos", Load(0.0, 0.0, -1.0, 0.0, 0.0), "D", dimensions, material)

end


function calculate_Mcrℓ_yy_pos(dimensions, material)

    return calculate_Rcr("Mcrℓ_yy_pos", Load(0.0, 0.0, -1.0, 0.0, 0.0), "L", dimensions, material)

end


function calculate_Mcrℓ_yy_neg(dimensions, material)

    return calculate_Rcr("Mcrℓ_yy_neg", Load(0.0, 0.0, 1.0, 0.0, 0.0), "L", dimensions, material)

end



# ---------------------------------------------------------------------------
# Net-section modes at a web hole
# ---------------------------------------------------------------------------

#Zero out the element thickness spanning a web hole (`hole_width`, centered at
#mid-depth) to model the net section in a finite strip buckling analysis (Moen &
#Schafer). `coordinates` is the gross (unperforated) section's centerline node
#geometry, `t` is the base thickness. Returns the net section coordinates and the
#per-element thickness vector: the nodes inside the hole span are replaced by two
#boundary nodes joined by ONE zero-thickness element, so the two parts of the section
#stay geometrically connected but carry no stiffness across the hole, and any node
#left within 0.1t of a boundary node is dropped so no near-degenerate strip is created
#(BucklingModeIdentification.zero_thickness_hole).
function calculate_net_section_element_thicknesses(coordinates, t, hole_width)

    X = coordinates.X
    Y = coordinates.Y

    web_lo, web_hi = CrossSectionGeometry.find_web_segment(X, Y)
    y_mid = (Y[web_lo] + Y[web_hi]) / 2

    Xc, Yc, t_elements = BMI.zero_thickness_hole(X, Y, t, y_mid - hole_width / 2, y_mid + hole_width / 2; tol = 0.1 * t)

    coordinates_net = (X=Xc, Y=Yc)

    return coordinates_net, t_elements

end


#Reduced web thickness model for distortional buckling at a web hole of length `Lh`
#(along the member axis), AISI S100-2024 Appendix 2 Section 2.3.3.3 (members with web
#holes; unchanged from S100-16 Eq. C-2.2.7-2, hence AISIS100.v16S3):
#tr = t*(1 - Lh/Lcrd)^(1/3), assigned to every web element (flanges/lips stay at the
#gross thickness). Returns the per-element thickness vector and tr.
function calculate_reduced_web_thicknesses(coordinates, t, Lh, Lcrd)

    Lh < Lcrd || error("Lh = $Lh must be less than Lcrd = $Lcrd for the AISI S100 Appendix 2 reduced thickness approach to apply")

    tr = AISIS100.v16S3.app2C2262(t, Lh, Lcrd)

    web_lo, web_hi = CrossSectionGeometry.find_web_segment(coordinates.X, coordinates.Y)
    i_start, i_end = minmax(web_lo, web_hi)

    t_elements = fill(t, length(coordinates.X) - 1)
    t_elements[i_start:(i_end - 1)] .= tr

    return t_elements, tr

end


#Elastic local buckling load, Pcrℓ, at the net section through a web hole of width
#`hole_width` (across the web, centered at mid-depth) and length `hole_length` (along
#the member axis): zero-thickness web element across the hole, Lcrℓ = min(hole_length,
#cFSM pure-local half-wavelength of the net section), conventional signature-curve
#value at that Lcrℓ (identify_local_buckling_at_hole).
function calculate_Pcrℓ_net_section(dimensions, material, hole_width, hole_length)

    label = "Pcrℓ net section"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    coordinates_net, t_elements = calculate_net_section_element_thicknesses(coordinates, dimensions.t, hole_width)

    results = identify_local_buckling_at_hole(coordinates_net.X, coordinates_net.Y, t_elements, load, material, hole_length; W = max_flat_width(coordinates.X, coordinates.Y))

    section = Section(label, material, dimensions, load, results)

    return section

end


#Elastic local buckling moment, Mcrℓ_xx, at the net section through a web hole of
#width `hole_width` and length `hole_length` (see calculate_Pcrℓ_net_section).
function calculate_Mcrℓ_xx_net_section(dimensions, material, hole_width, hole_length)

    label = "Mcrℓ_xx net section"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    coordinates_net, t_elements = calculate_net_section_element_thicknesses(coordinates, dimensions.t, hole_width)

    results = identify_local_buckling_at_hole(coordinates_net.X, coordinates_net.Y, t_elements, load, material, hole_length; W = max_flat_width(coordinates.X, coordinates.Y))

    section = Section(label, material, dimensions, load, results)

    return section

end


#Elastic distortional buckling load, Pcrd, at a section with a web hole of length
#`Lh` (measured along the member axis, not across the web - contrast with the
#`hole_width` argument of the local buckling net-section functions, which is a
#web-depth dimension), using the reduced web thickness approach of AISI S100-2024
#Appendix 2 Section 2.3.3.3 (see calculate_reduced_web_thicknesses). Lcrd is the gross
#section's identified distortional buckling half-wavelength (calculate_Pcrd), and the
#reduced-thickness model is evaluated at that same Lcrd (method = :specified_length).
#
#This numerically approximates the closed-form provision rather than reproducing it
#exactly: the provision derives tr specifically so tr^3/t^3 scales kϕwe (the web's
#elastic distortional stiffness term, cubic in t), leaving kϕwg (web geometric
#stiffness, linear in t) and the flange terms at full thickness. Assigning tr to the
#web elements in a finite strip model scales both the web's elastic AND geometric
#stiffness contributions, not kϕwe alone.
function calculate_Pcrd_net_section(dimensions, material, Lh)

    label = "Pcrd net section (reduced web thickness)"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    gross = calculate_Pcrd(dimensions, material)
    Lcrd = gross.results.Lcr

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    t_elements, tr = calculate_reduced_web_thicknesses(coordinates, dimensions.t, Lh, Lcrd)

    results = evaluate_at_length(coordinates.X, coordinates.Y, t_elements, load, material, Lcrd; mode = :distortional, Lcr_estimate = gross.results.identification.Lcr_estimate)

    section = Section(label, material, dimensions, load, results)

    return section

end


#Elastic distortional buckling moment, Mcrd_xx, at a section with a web hole of
#length `Lh` (along the member axis), same reduced web thickness approach as
#calculate_Pcrd_net_section with Lcrd the gross section's identified distortional
#buckling half-wavelength UNDER BENDING (calculate_Mcrd_xx).
function calculate_Mcrd_xx_net_section(dimensions, material, Lh)

    label = "Mcrd_xx net section (reduced web thickness)"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    gross = calculate_Mcrd_xx(dimensions, material)
    Lcrd = gross.results.Lcr

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    t_elements, tr = calculate_reduced_web_thicknesses(coordinates, dimensions.t, Lh, Lcrd)

    results = evaluate_at_length(coordinates.X, coordinates.Y, t_elements, load, material, Lcrd; mode = :distortional, Lcr_estimate = gross.results.identification.Lcr_estimate)

    section = Section(label, material, dimensions, load, results)

    return section

end



# ---------------------------------------------------------------------------
# Collections
# ---------------------------------------------------------------------------

function calculate_properties(dimensions, material)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates
    t = dimensions.t
    section_properties = calculate_section_properties(coordinates, t)

    mode_labels = ["Pcrℓ", "Pcrd", "Mcrd_xx", "Mcrℓ_xx", "Mcrd_yy_pos",  "Mcrℓ_yy_pos",  "Mcrℓ_yy_neg"]

    mode_loads = [Load(1.0, 0.0, 0.0, 0.0, 0.0),
                    Load(1.0, 0.0, 0.0, 0.0, 0.0),
                    Load(0.0, 1.0, 0.0, 0.0, 0.0),
                    Load(0.0, 1.0, 0.0, 0.0, 0.0),
                    Load(0.0, 0.0, -1.0, 0.0, 0.0),
                    Load(0.0, 0.0, -1.0, 0.0, 0.0),
                    Load(0.0, 0.0, 1.0, 0.0, 0.0)]


    mode_types = ["L", "D", "D", "L", "D", "L", "L"]

    values = [section_properties; [calculate_Rcr(mode_labels[i], mode_loads[i], mode_types[i], dimensions, material) for i in eachindex(mode_labels)]]  # one value per label, same order as mode_labels

    fields = ["section"; mode_labels]

    properties = NamedTuple{Tuple(Symbol.(fields))}(Tuple(values))

    return properties

end


function calculate_all_buckling_modes(dimensions, material)


    results = (Pcrℓ = calculate_Pcrℓ(dimensions, material),
    Pcrd = calculate_Pcrd(dimensions, material),
    Mcrd_xx = calculate_Mcrd_xx(dimensions, material),
    Mcrℓ_xx = calculate_Mcrℓ_xx(dimensions, material),
    Mcrd_yy_pos = calculate_Mcrd_yy_pos(dimensions, material),
    Mcrℓ_yy_pos = calculate_Mcrℓ_yy_pos(dimensions, material),
    Mcrℓ_yy_neg = calculate_Mcrℓ_yy_neg(dimensions, material))


end


#One-line description of how a Results' critical point was located.
function describe(results::Results)

    id = results.identification
    id === nothing && return "minimum of the given half-wavelengths: Rcr = $(results.Rcr) at Lcr = $(results.Lcr)"

    what = id.method == :signature_trough ? "first trough of the conventional signature curve" :
           id.method == :cFSM_pure_mode ? "no trough in the search range: conventional value at the cFSM pure-$(id.cFSM_subspace) half-wavelength" :
           "conventional value at a specified half-wavelength"

    range = id.W === nothing ? "search range $(round.(id.search_range, digits = 3))" : "W = $(round(id.W, digits = 3)), search range $(round.(id.search_range, digits = 3))"
    cfsm = id.Lcr_cFSM === nothing ? "" : "; cFSM pure-$(id.cFSM_subspace) minimum $(round(id.Rcr_cFSM, digits = 4)) at L = $(round(id.Lcr_cFSM, digits = 4))"
    p = id.participation

    return "$(id.mode) buckling, $(what): Rcr = $(round(results.Rcr, digits = 4)) at Lcr = $(round(results.Lcr, digits = 4)) ($(range)$(cfsm); participation G $(round(p.G, digits = 1)) % D $(round(p.D, digits = 1)) % L $(round(p.L, digits = 1)) % O $(round(p.O, digits = 1)) %)"

end





end # module CeeSectionBuckling
