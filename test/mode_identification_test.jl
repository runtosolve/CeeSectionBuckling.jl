# Buckling mode identification on the actual rounded-corner geometry
# (BucklingModeIdentification.jl): gross and net-section modes of the MRI 362S400-68
# sample-calculation stud (S400 flange: B = 4.0 in., lip 1.0 in.; 68 mil), which has a
# clear local trough, and a check that every identified point is what the procedure
# says it is. Run from an environment that has CeeSectionBuckling (e.g.
# ~/.julia/dev/SFIASections/generate):
#   julia +1.12 --startup-file=no --project=. mode_identification_test.jl
using CeeSectionBuckling, Test
const CSB = CeeSectionBuckling
const CUFSM = CeeSectionBuckling.CUFSM

E = 29500.0
ν = 0.30

#362S400-68
t = 0.0713
L = 1.0
B = 4.0
H = 3.62
r = 0.1070

material = CSB.Material(E, ν)
dimensions = CSB.Dimensions(t, L, B, H, r)

hole_width = 1.5   #in., SFIA punchout width for H > 2.5 in.
hole_length = 4.0  #in.

function check_identified(section)
    res = section.results
    id = res.identification
    @test id !== nothing
    @test res.model.lengths == [res.Lcr]
    @test CUFSM.Tools.get_load_factor(res.model, 1)[1] ≈ res.Rcr
    @test id.method in (:signature_trough, :cFSM_pure_mode, :hole_length, :specified_length)
    if id.method == :signature_trough
        i = findfirst(==(res.Lcr), id.lengths)
        @test i !== nothing && i in id.troughs
        @test id.curve[i] ≈ res.Rcr rtol = 1e-6
        p = id.participation
        @test CSB.dominant(p) == (id.mode == :local ? :L : :D)       #the trough is of the mode it is reported as
        if id.mode == :local
            @test i == first(id.troughs)                              #FIRST L-dominated trough in the range
            @test id.search_range[1] <= res.Lcr <= id.search_range[2]
        end
    elseif id.method == :cFSM_pure_mode
        @test res.Lcr == id.Lcr_cFSM
        @test res.Rcr ≈ id.Rcr_at_Lcr_cFSM
        @test res.Rcr <= id.Rcr_cFSM * (1 + 1e-6)                   #constrained minimum is an upper bound
    end
    if id.Lcr_cFSM !== nothing && id.method != :specified_length
        @test id.Rcr_at_Lcr_cFSM <= id.Rcr_cFSM * (1 + 1e-6)
    end
    p = id.participation
    @test p.G + p.D + p.L + p.O ≈ 100 atol = 1e-6
    println(rpad(section.label, 42), CSB.describe(res))
    return res
end

@testset "362S400-68 gross section" begin
    Pcrℓ = check_identified(CSB.calculate_Pcrℓ(dimensions, material))
    Mcrℓ = check_identified(CSB.calculate_Mcrℓ_xx(dimensions, material))
    Pcrd = check_identified(CSB.calculate_Pcrd(dimensions, material))
    Mcrd = check_identified(CSB.calculate_Mcrd_xx(dimensions, material))
    @test Pcrℓ.identification.mode == :local && Mcrℓ.identification.mode == :local
    @test Pcrd.identification.mode == :distortional && Mcrd.identification.mode == :distortional
    @test Mcrℓ.identification.participation.L > 50           #the local trough is L-dominated
    @test Mcrd.identification.participation.D > 50           #the distortional trough is D-dominated
    @test Mcrℓ.identification.W ≈ B - 2 * (r + t) atol = 0.02 #S400: the flange flat is the widest element
    #MRI report Revision 2 (September 2026) values for this section, 9-point search: Mcrℓ = 55.88 kip-in.
    @test Mcrℓ.Rcr ≈ 55.88 rtol = 0.01
    #Lcrd is the trough of the gross signature curve (about 43 in.), not the App. 2 estimate (29.3 in.)
    @test Mcrd.identification.Lcr_estimate ≈ 29.33 rtol = 0.01
    @test Mcrd.Lcr > 1.2 * Mcrd.identification.Lcr_estimate
    @test Mcrd.Rcr < 75.79                                   #below the value at the App. 2 estimate
end

@testset "362S400-68 net section (1.5 x 4.0 in. punchout)" begin
    Mcrℓ_net = check_identified(CSB.calculate_Mcrℓ_xx_net_section(dimensions, material, hole_width, hole_length))
    Pcrℓ_net = check_identified(CSB.calculate_Pcrℓ_net_section(dimensions, material, hole_width, hole_length))
    Mcrd_net = check_identified(CSB.calculate_Mcrd_xx_net_section(dimensions, material, hole_length))
    Pcrd_net = check_identified(CSB.calculate_Pcrd_net_section(dimensions, material, hole_length))
    @test count(==(0.0), Mcrℓ_net.model.elem[:, 4]) == 1     #one zero-thickness element across the hole
    @test Mcrℓ_net.Rcr ≈ 51.22 rtol = 0.01                   #MRI report Revision 2 value 51.23
    #net-section local half-wavelength is min(Lh, Lcrℓ,cFSM): here the cFSM minimum (3.89 in.) is shorter than the hole
    @test Mcrℓ_net.identification.method == :cFSM_pure_mode
    @test Mcrℓ_net.Lcr ≈ 3.89 rtol = 0.01
    @test Mcrℓ_net.Lcr <= hole_length && Pcrℓ_net.Lcr <= hole_length
    cg = CSB.get_section_coordinates(dimensions)
    Wg = CSB.max_flat_width(cg.X, cg.Y)
    @test Mcrℓ_net.identification.W ≈ Wg                     #net models are searched over the gross section's [0.5W, 2W]
    @test Mcrℓ_net.identification.search_range == (0.5 * Wg, 2.0 * Wg)
    @test Mcrd_net.identification.method == :specified_length
    @test Mcrd_net.Lcr == CSB.calculate_Mcrd_xx(dimensions, material).results.Lcr
    @test any(<(t), Mcrd_net.model.elem[:, 4]) && all(>(0.0), Mcrd_net.model.elem[:, 4])   #reduced web thickness
end

@testset "deep web net section: local half-wavelength limited to the hole length" begin
    #1200S162-33 with the 1.5 x 4.0 in. punchout: the pure-local cFSM minimum of the net
    #section is longer than the 4.0 in. hole, so Lcrℓ = Lh and the conventional value there is reported
    d = CSB.Dimensions(0.0346, 0.5, 1.625, 12.0, 0.0765)
    net = check_identified(CSB.calculate_Mcrℓ_xx_net_section(d, material, hole_width, hole_length))
    @test net.identification.method == :hole_length
    @test net.Lcr == hole_length
    @test net.identification.Lcr_cFSM > hole_length
    @test net.Rcr ≈ 14.50 rtol = 0.01
end

@testset "identification helpers" begin
    c = CSB.get_section_coordinates(dimensions)
    @test CSB.max_flat_width(c.X, c.Y) ≈ B - 2 * (r + t) atol = 0.02
    @test CSB.troughs([3.0, 2.0, 1.0, 2.0, 1.5, 1.0, 0.5]) == [3]
    @test isempty(CSB.troughs([3.0, 2.0, 1.0]))
    #low-level sweeps keep their minimum-of-the-given-lengths behaviour, identification = nothing
    res = CSB.calculate_buckling_properties(c, dimensions, CSB.Load(0.0, 1.0, 0.0, 0.0, 0.0), material, [3.0, 3.5, 4.0])
    @test res.identification === nothing && res.Lcr in (3.0, 3.5, 4.0)
end

@testset "deep web, narrow thick flange: distortional trough inside [0.5W, 2W]" begin
    #1150S137-43: the first trough of the bending signature curve inside [0.5W, 2W] is
    #distortional-dominated (D ~ 49 %, L ~ 49 %) and the window about the Appendix 2 Lcrd
    #estimate contains the web's local trough; the participation must keep them apart
    d = CSB.Dimensions(0.0451, 0.375, 1.375, 11.5, 0.0712)
    Mcrℓ = check_identified(CSB.calculate_Mcrℓ_xx(d, material))
    Mcrd = check_identified(CSB.calculate_Mcrd_xx(d, material))
    @test CSB.dominant(Mcrℓ.identification.participation) == :L
    @test CSB.dominant(Mcrd.identification.participation) == :D
    @test Mcrd.Lcr > Mcrℓ.Lcr                                   #distortional half-wavelength beyond the local one
    @test !(Mcrℓ.Rcr ≈ Mcrd.Rcr)                                #not the same trough read twice
    Pcrd = check_identified(CSB.calculate_Pcrd(d, material))
    @test CSB.dominant(Pcrd.identification.participation) == :D
end

@testset "cFSM pure-local fallback (search range with no trough)" begin
    c = CSB.get_section_coordinates(dimensions)
    load = CSB.Load(0.0, 1.0, 0.0, 0.0, 0.0)
    #[0.5W, 0.7W] lies on the falling branch of the curve (the trough is near 0.95W)
    res = CSB.identify_local_buckling(c.X, c.Y, t, load, material; range_factors = (0.5, 0.7), n_range = 9)
    @test isempty(res.identification.troughs)
    @test res.identification.method == :cFSM_pure_mode
    @test res.Lcr == res.identification.Lcr_cFSM
    @test res.Rcr ≈ res.identification.Rcr_at_Lcr_cFSM
    @test res.Rcr <= res.identification.Rcr_cFSM * (1 + 1e-6)
    @test res.model.lengths == [res.Lcr]
    println(rpad("forced fallback", 42), CSB.describe(res))
    #without the constrained analysis the same range is an error
    @test_throws ErrorException CSB.identify_local_buckling(c.X, c.Y, t, load, material; range_factors = (0.5, 0.7), n_range = 9, cfsm = false)
end
