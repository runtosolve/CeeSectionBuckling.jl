using CeeSectionBuckling

E = 29500.0
ν = 0.30 

t = 0.102 
L = 0.91
B = 2.5
D = 6.0
r = 0.1875


material = CeeSectionBuckling.Material(E, ν)
dimensions = CeeSectionBuckling.Dimensions(t, L, B, D, r)

section = CeeSectionBuckling.calculate_Pcrd(dimensions, material)