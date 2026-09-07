using CeeSectionBuckling, ShowCUFSM, CairoMakie

E = 29500.0
ν = 0.30

#362S162-33
t = 0.0346
L = 0.5
B = 1.625
H = 3.625
r = 0.0765

material = CeeSectionBuckling.Material(E, ν)
dimensions = CeeSectionBuckling.Dimensions(t, L, B, H, r)

Lh = 2.0  #hole length along the member axis, in.

gross_section = CeeSectionBuckling.calculate_Pcrd(dimensions, material)
net_section = CeeSectionBuckling.calculate_Pcrd_net_section(dimensions, material, Lh)

println("Pcrd, gross section = ", gross_section.results.Rcr, " kips at Lcr = ", gross_section.results.Lcr, " in.")
println("Pcrd, net section (", Lh, " in. long web hole) = ", net_section.results.Rcr, " kips at Lcr = ", net_section.results.Lcr, " in.")

#distortional buckling mode shape at the net section, via
#ShowCUFSM.cross_section_mode_shape_info (picks the governing length automatically;
#plotted per element below, with line width proportional to each element's actual
#thickness, so the web's reduced thickness tr shows up as a visibly thinner run than
#the flanges/lips).
#
#ShowCUFSM.minimum_mode_shape wraps this same call with a convenience line-merging
#step (LinesCurvesNodes.find_linesegments/combine_points_into_linesegments) that can
#throw a BoundsError on certain thickness-run patterns, so this test calls the
#lower-level, per-element routine directly instead (see Pcrl_net_section_test.jl).
model = net_section.results.model
t_elements = model.elem[:, 4]
eig = 1
deformation_scale = [1.0, 1.0]

Pcr = ShowCUFSM.get_load_factor(model, eig)
mode_index = argmin(Pcr)
mode = model.shapes[mode_index][:, eig]

n = fill(10, length(t_elements))
cross_section_coords, Δ, _ = ShowCUFSM.cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)

f = Figure(size = (700, 800))
ax = Axis(f[1, 1], aspect = DataAspect(),
    title = "362S162-33 net section distortional buckling, $(Lh) in. long web hole\nPcrd,net = $(round(net_section.results.Rcr, digits = 2)) kips, Lcr = $(round(net_section.results.Lcr, digits = 2)) in.")
hidedecorations!(ax)
hidespines!(ax)

thickness_scale = 300.0

for i in eachindex(cross_section_coords)

    X = [p[1] for p in cross_section_coords[i]]
    Y = [p[2] for p in cross_section_coords[i]]
    ΔX = [d[1] for d in Δ[i]]
    ΔY = [d[2] for d in Δ[i]]

    linewidth = max(t_elements[i] * thickness_scale, 1.0)

    lines!(ax, X, Y, color = (:gray, 0.5), linestyle = :dash, linewidth = 1)
    lines!(ax, X .+ deformation_scale[1] .* ΔX, Y .+ deformation_scale[2] .* ΔY, color = :blue, linewidth = linewidth)

end

save(joinpath(@__DIR__, "Pcrd_net_section_mode_shape.png"), f)

f
