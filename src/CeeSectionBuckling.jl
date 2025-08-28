module CeeSectionBuckling


using CrossSectionGeometry, CUFSM, AISIS100 


struct Material

    E
    ν

end

struct Dimensions

    t 

    L 
    B 
    D 
    r  #inside radius 

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
    D, 
    r  #inside radius 
    ) = dimensions 

    section_dimensions = [L, B, D, B, L]
    r = [r+t, r+t, r+t, r+t]
    n = [3, 3, 3, 3, 3]
    n_r = [3, 3, 3, 3];
    θ = [π/2, π, -π/2, 0.0, π/2]

    coordinates = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, r, n_r, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3] - section_dimensions[1]))

    X = [coordinates.centerline_node_XY[i][1] for i in eachindex(coordinates.centerline_node_XY)]
    Y = [coordinates.centerline_node_XY[i][2] for i in eachindex(coordinates.centerline_node_XY)]

    coordinates = (X=X, Y=Y)

    return coordinates

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
    D, 
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
    ho = D
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



function calculate_Pcrℓ(dimensions, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    if B > D
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * D, 1.5 * D, 9)
    end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Pcrℓ(dimensions, coordinates, material)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    if B > D
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.5 * D, 1.5 * D, 9)
    end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Pcrℓ(dimensions, coordinates, material, lengths)

    label = "Pcrℓ"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    # B = dimensions.B 
    # D = dimensions.D

    # if B > D
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.5 * D, 1.5 * D, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

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



function calculate_Pcrd(dimensions, coordinates, material, lengths)

    label = "Pcrd"

    load = Load(1.0, 0.0, 0.0, 0.0, 0.0)

    # load_type = "P"
    # lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

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


function calculate_Mcrd_xx(dimensions, coordinates, material)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrd_xx(dimensions, coordinates, material, lengths)

    label = "Mcrd_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    load_type = "M"
    # lengths = [calculate_Lcrd(dimensions, material, load_type)]

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Mcrℓ_xx(dimensions, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    if B > D / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * D, 1.25 * D, 9)
    end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end



function calculate_Mcrℓ_xx(dimensions, coordinates, material)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    if B > D / 2
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * D, 1.25 * D, 9)
    end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_xx(dimensions, coordinates, material, lengths)

    label = "Mcrℓ_xx"

    load = Load(0.0, 1.0, 0.0, 0.0, 0.0)

    # B = dimensions.B 
    # D = dimensions.D

    # if B > D / 2
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

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

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrd_yy_pos(dimensions, coordinates, material)

    label = "Mcrd_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    load_type = "M"
    Lcrd = calculate_Lcrd(dimensions, material, load_type)
    
    lengths = range(0.5 * Lcrd, 1.5 * Lcrd, 9)

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrℓ_yy_pos(dimensions, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    # lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)

    if B > 2 * D 
        lengths = range(0.5 * B, 1.5 * B, 9)
    else
        lengths = range(0.25 * D, 1.25 * D, 9)
    end

    # B = dimensions.B 
    # D = dimensions.D

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrℓ_yy_pos(dimensions, coordinates, material)

    label = "Mcrℓ_yy_pos"

    load = Load(0.0, 0.0, -1.0, 0.0, 0.0)

    B = dimensions.B 
    # D = dimensions.D
    L = dimensions.L

    # lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)

    lengths = range(L, 1.2 * B, 9)

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end


    # B = dimensions.B 
    # D = dimensions.D

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end






function calculate_Mcrℓ_yy_neg(dimensions, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D


    lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)



    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




function calculate_Mcrℓ_yy_neg(dimensions, coordinates, material)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)



    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end


function calculate_Mcrℓ_yy_neg(dimensions, coordinates, material, lengths)

    label = "Mcrℓ_yy_neg"

    load = Load(0.0, 0.0, 1.0, 0.0, 0.0)

    B = dimensions.B 
    D = dimensions.D

    # if B > 2 * D 
    #     lengths = range(0.5 * B, 1.5 * B, 9)
    # else
    #     lengths = range(0.25 * D, 1.25 * D, 9)
    # end

    # lengths = range(0.5*maximum([B, D]), 2.0*maximum([B, D]), 7)



    # coordinates = get_section_coordinates(dimensions)

    results = calculate_buckling_properties(coordinates, dimensions, load, material, lengths)

    section = Section(label, material, dimensions, load, results)

    return section 

end




end # module CeeSectionBuckling
