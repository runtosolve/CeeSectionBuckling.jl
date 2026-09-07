module CeeSectionBuckling


using CrossSectionGeometry, CUFSM, AISIS100, SectionProperties, OpenSectionBuckling


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

struct Results 

    model 
    Lcr
    Rcr 

end


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


function calculate_constrained_buckling_properties(straight_corner_coordinates, dimensions, loads, material, mode_type)

    (;E,
    ν
    ) = material 

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    t = dimensions.t



    ###

    # flat_mesh_size_goal = 0.5
    # corner_mesh_size_goal = π/6
    # centerline_radius = r + t/2


    #coordinates coming in are straight corner X, Y

    coordinates = (straight=straight_corner_coordinates, rounded=straight_corner_coordinates) #rounded not used for now 

    material = (E=E, ν=ν)

    loads = (P=P, Mxx=Mxx, Mzz=Mzz, M11=M11, M22=M22)

    # mode_type = ["D"]
    model = OpenSectionBuckling.properties(coordinates, material, loads, t, mode_type)

    return model 

end


function calculate_section_properties(coordinates, t)

    X = coordinates.X 
    Y = coordinates.Y

    t_all = fill(t, length(X)-1)

    section_properties = SectionProperties.open_thin_walled(X, Y, t_all)

    return section_properties 

end








function calculate_buckling_properties(coordinates, dimensions, loads, material, lengths)

    (;E,
    ν
    ) = material 

    (;P,
    Mxx,
    Mzz,
    M11,
    M22,
    ) = loads

    t = dimensions.t

    constraints = []
    springs = []
    
    neigs = 1

    model = CUFSM.Tools.open_section_analysis(coordinates.X, coordinates.Y, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)
    Rcr = minimum(Rcr_curve) 
    
    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results 

end

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


function calculate_Rcr(label, load, mode_type, dimensions, material)

    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    lengths = range(0.75 * Lcr_cFSM, 1.25 * Lcr_cFSM, 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Pcrℓ(dimensions, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    H = dimensions.H

    if B > H
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * H, 1.5 * H, 9)
    end

    coordinates = get_section_coordinates(dimensions) #centerline coordinates 

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    #catch indistinct local buckling minimum here 
    mode_type = ["L"]
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    # println("Lcr_cFSM", Lcr_cFSM)
    if !isapprox(results.Lcr, Lcr_cFSM, rtol=0.05)

        # println("results.Lcr", results.Lcr)
        # println("Lcr_cFSM", Lcr_cFSM)
        lengths = [Lcr_cFSM]
        results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

        section = Section(label, material, dimensions, load, results)

    else

        section = Section(label, material, dimensions, load, results)

    end
    

    return section 

end




#Zero out the element thicknesses spanning a web hole (`hole_width`, centered at
#mid-depth) to model the net section in a finite strip buckling analysis.
#`coordinates` is the gross (unperforated) section's centerline node geometry, `t` is
#the base thickness. Returns the augmented net section coordinates and the
#per-element thickness vector (zero across the hole span).
function calculate_net_section_element_thicknesses(coordinates, t, hole_width)

    X = coordinates.X
    Y = coordinates.Y

    Xc, Yc, idx_lo, idx_hi = CrossSectionGeometry.insert_hole_span_boundary_nodes(X, Y, hole_width)

    #Any original mesh nodes strictly between the hole boundary nodes are interior to
    #the hole once their flanking elements are zeroed out, leaving them connected only
    #to zero-thickness elements on both sides. With no real material tying them to the
    #rest of the model, they're floating, unconstrained degrees of freedom that
    #corrupt the eigenvalue analysis (a spurious near-zero-energy "mode" of just those
    #nodes swamps the real buckling mode). They're colinear with idx_lo/idx_hi (all on
    #the same flat), so drop them and span the hole with a single zero-thickness
    #element.
    keep = vcat(1:idx_lo, idx_hi:length(Xc))
    Xc = Xc[keep]
    Yc = Yc[keep]
    idx_hi = idx_lo + 1

    t_elements = fill(t, length(Xc) - 1)
    t_elements[idx_lo] = 0.0

    coordinates_net = (X=Xc, Y=Yc)

    return coordinates_net, t_elements

end


#Same as `calculate_buckling_properties`, but takes the element thicknesses directly
#(scalar or per-element vector) instead of pulling a single thickness from
#`dimensions`, so it can be used with the net section's zeroed-out element
#thicknesses at a web hole.
function calculate_net_section_buckling_properties(coordinates, t, loads, material, lengths)

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

    model = CUFSM.Tools.open_section_analysis(coordinates.X, coordinates.Y, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    eig = 1
    Rcr_curve = CUFSM.Tools.get_load_factor(model, eig)
    Rcr = minimum(Rcr_curve)

    index = argmin(Rcr_curve)
    Lcr = lengths[index]

    results = Results(model, Lcr, Rcr)

    return results

end


#Elastic local buckling load, Pcrℓ, at the net section through a web hole of width
#`hole_width` centered at mid-depth. The web elements spanning the hole are modeled
#with zero thickness in the finite strip analysis, per calculate_net_section_element_thicknesses.
function calculate_Pcrℓ_net_section(dimensions, material, hole_width)

    label = "Pcrℓ net section"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    #A hole can localize buckling over a shorter half-wavelength than the gross
    #section's flat-width-controlled local mode, so search down to a fraction of the
    #hole width as well.
    lengths = range(0.25 * minimum([B, H, hole_width]), 1.5 * maximum([B, H]), 9)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    coordinates_net, t_elements = calculate_net_section_element_thicknesses(coordinates, dimensions.t, hole_width)

    results = calculate_net_section_buckling_properties(coordinates_net, t_elements, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end


function calculate_Pcrd(dimensions, material)

    label = "Pcrd"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    load_type = "P"
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end


#Elastic distortional buckling load, Pcrd, at a section with a web hole of length
#`Lh` (measured along the member axis, not across the web - contrast with the
#`hole_width` argument of calculate_Pcrℓ_net_section/calculate_Mcrℓ_xx_net_section,
#which is a web-depth dimension for their local buckling finite strip approach), using
#the reduced web thickness approach of AISI S100-2024 Appendix 2 (2.3.3.6, "Members
#with a Web Hole" - carried over unchanged from S100-16, hence AISIS100.v16S3):
#tr = t*(1 - Lh/Lcrd)^(1/3), where Lcrd is the gross section's distortional buckling
#half-wavelength. `tr` is assigned to every web element (flanges/lips stay at the
#gross thickness) over the section's full length, and the finite strip eigenvalue
#analysis is run at Lcrd, same as calculate_Pcrd's gross-section calculation.
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

    load_type = "P"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)

    Lh < Lcrd || error("Lh = $Lh must be less than Lcrd = $Lcrd for the AISI S100 Appendix 2 reduced thickness approach to apply")

    tr = AISIS100.v16S3.app2C2262(dimensions.t, Lh, Lcrd)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    web_lo, web_hi = CrossSectionGeometry.find_web_segment(coordinates.X, coordinates.Y)
    i_start, i_end = minmax(web_lo, web_hi)

    t_elements = fill(dimensions.t, length(coordinates.X) - 1)
    t_elements[i_start:(i_end - 1)] .= tr

    lengths = [Lcrd]

    results = calculate_net_section_buckling_properties(coordinates, t_elements, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end




function calculate_Mcrd_xx(dimensions, material)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end


#Elastic distortional buckling moment, Mcrd_xx, at a section with a web hole of
#length `Lh` (measured along the member axis - see the comment on
#calculate_Pcrd_net_section for the contrast with the `hole_width` used by the local
#buckling net-section functions), using the same AISI S100-2024 Appendix 2 reduced
#web thickness approach: tr = t*(1 - Lh/Lcrd)^(1/3), with Lcrd here the gross
#section's distortional buckling half-wavelength UNDER BENDING (load_type = "M").
#`tr` is assigned to every web element (flanges/lips stay at the gross thickness)
#and the finite strip eigenvalue analysis is run at Lcrd, same as calculate_Mcrd_xx's
#gross-section calculation. See calculate_Pcrd_net_section for the same caveat about
#this being a practical finite-strip approximation of the closed-form provision
#(which only substitutes tr into kϕwe), not an exact reproduction of it.
function calculate_Mcrd_xx_net_section(dimensions, material, Lh)

    label = "Mcrd_xx net section (reduced web thickness)"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)

    Lh < Lcrd || error("Lh = $Lh must be less than Lcrd = $Lcrd for the AISI S100 Appendix 2 reduced thickness approach to apply")

    tr = AISIS100.v16S3.app2C2262(dimensions.t, Lh, Lcrd)

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    web_lo, web_hi = CrossSectionGeometry.find_web_segment(coordinates.X, coordinates.Y)
    i_start, i_end = minmax(web_lo, web_hi)

    t_elements = fill(dimensions.t, length(coordinates.X) - 1)
    t_elements[i_start:(i_end - 1)] .= tr

    lengths = [Lcrd]

    results = calculate_net_section_buckling_properties(coordinates, t_elements, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end



function calculate_Mcrℓ_xx(dimensions, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    #catch indistinct local buckling minimum here 
    mode_type = ["L"]
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    if !isapprox(results.Lcr, Lcr_cFSM, rtol=0.05)

        # println("results.Lcr", results.Lcr)
        # println("Lcr_cFSM", Lcr_cFSM)
        lengths = [Lcr_cFSM]
        results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

        section = Section(label, material, dimensions, load, results)

    else

        section = Section(label, material, dimensions, load, results)

    end

    return section

end


#Elastic local buckling moment, Mcrℓ_xx, at the net section through a web hole of
#width `hole_width` centered at mid-depth. The web elements spanning the hole are
#modeled with zero thickness in the finite strip analysis, per
#calculate_net_section_element_thicknesses. As with calculate_Pcrℓ_net_section, this
#skips the indistinct-local-buckling-minimum cFSM correction used for the gross
#section, since that constrained modal classification doesn't have a clean analogue
#once part of the web is zero-thickness.
function calculate_Mcrℓ_xx_net_section(dimensions, material, hole_width)

    label = "Mcrℓ_xx net section"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B
    H = dimensions.H

    if B > H / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    coordinates = get_section_coordinates(dimensions) #centerline coordinates, gross section

    coordinates_net, t_elements = calculate_net_section_element_thicknesses(coordinates, dimensions.t, hole_width)

    results = calculate_net_section_buckling_properties(coordinates_net, t_elements, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section

end




#this is not the traditional strong axis distortional buckling, still important though in some cases for weak axis bending
function calculate_Mcrd_yy_pos(dimensions, material)

    label = "Mcrd_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    load_type = "M"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(0.5 * Lcrd, 1.5 * Lcrd, 9)

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    # section = Section(label, material, dimensions, load, results)

    #catch indistinct local buckling minimum here 
    mode_type = ["D"]
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    # println("Lcr_cFSM", Lcr_cFSM)
    if !isapprox(results.Lcr, Lcr_cFSM, rtol=0.05)

        # println("results.Lcr", results.Lcr)
        # println("Lcr_cFSM", Lcr_cFSM)
        lengths = [Lcr_cFSM]
        results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

        section = Section(label, material, dimensions, load, results)

    else

        section = Section(label, material, dimensions, load, results)

    end


    return section 

end



function calculate_Mcrℓ_yy_pos(dimensions, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B 
    H = dimensions.H

    # lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)

    if B > 2 * H 
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * H, 1.25 * H, 9)
    end

    # B = dimensions.B 
    # H = dimensions.H

    # if B > 2 * H 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    #catch indistinct local buckling minimum here 
    mode_type = ["L"]
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    if !isapprox(results.Lcr, Lcr_cFSM, rtol=0.05)

        # println("results.Lcr", results.Lcr)
        # println("Lcr_cFSM", Lcr_cFSM)
        lengths = [Lcr_cFSM]
        results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

        section = Section(label, material, dimensions, load, results)

    else

        section = Section(label, material, dimensions, load, results)

    end

    return section 

end






function calculate_Mcrℓ_yy_neg(dimensions, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B 
    H = dimensions.H


    lengths = range(0.5*maximum([B, H]), 2.0*maximum([B, H]), 7)



    # if B > 2 * H 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * H, 1.25 * H, 9)
    # end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    # section = Section(label, material, dimensions, load, results)

    #catch indistinct local buckling minimum here 
    mode_type = ["L"]
    coordinates = get_straight_corner_section_coordinates(dimensions) #straight corner coordinates 
    model = calculate_constrained_buckling_properties(coordinates, dimensions, load, material, mode_type)

    Lcr_cFSM = model.curve[1][1]

    # println("Lcr_cFSM", Lcr_cFSM)
    if !isapprox(results.Lcr, Lcr_cFSM, rtol=0.05)

        # println("results.Lcr", results.Lcr)
        # println("Lcr_cFSM", Lcr_cFSM)
        lengths = [Lcr_cFSM]
        results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

        section = Section(label, material, dimensions, load, results)

    else

        section = Section(label, material, dimensions, load, results)

    end


    return section 

end




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





end # module CeeSectionBuckling
