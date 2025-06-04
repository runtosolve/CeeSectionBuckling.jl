using CeeSectionBuckling


E = 29500.0
ν = 0.30 

t = 0.060
L = 3/8
B = 1.0
D = 2.0
r = 2 * t


material = CeeSectionBuckling.Material(E, ν)
dimensions = CeeSectionBuckling.Dimensions(t, L, B, D, r)


results = CeeSectionBuckling.calculate_Mcrℓ_yy_neg(dimensions, material)