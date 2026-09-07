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

hole_width = 1.5

gross_section = CeeSectionBuckling.calculate_Pcrℓ(dimensions, material)
net_section = CeeSectionBuckling.calculate_Pcrℓ_net_section(dimensions, material, hole_width)

println("Pcrℓ, gross section = ", gross_section.results.Rcr, " kips at Lcr = ", gross_section.results.Lcr, " in.")
println("Pcrℓ, net section (", hole_width, " in. hole) = ", net_section.results.Rcr, " kips at Lcr = ", net_section.results.Lcr, " in.")

#local buckling mode shape at the net section, via ShowCUFSM.cross_section_mode_shape_info
#(picks the governing length automatically; plotted per element below, with line
#width proportional to each element's actual thickness, so the zero-thickness hole
#element shows up as a hairline gap in the buckled shape).
#
#ShowCUFSM.minimum_mode_shape wraps this same call with a convenience line-merging
#step (LinesCurvesNodes.find_linesegments/combine_points_into_linesegments) that
#throws a BoundsError whenever a thickness run is exactly one raw element long, which
#is exactly what a hole spanning a single mesh element produces here - so this test
#calls the lower-level, per-element routine directly instead.
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
    title = "362S162-33 net section local buckling, $(hole_width) in. web hole\nPcrℓ,net = $(round(net_section.results.Rcr, digits = 2)) kips, Lcr = $(round(net_section.results.Lcr, digits = 2)) in.")
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

save(joinpath(@__DIR__, "Pcrl_net_section_mode_shape.png"), f)

f
