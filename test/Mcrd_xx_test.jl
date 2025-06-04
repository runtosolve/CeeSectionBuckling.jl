using CeeSectionBuckling

E = 29500.0
ν = 0.30 

t = 0.102 
L = 0.91
B = 2.5
D = 6.0
r = 0.1875

input = CeeSectionBuckling.SectionInput(E, ν, t, L, B, D, r)
output = CeeSectionBuckling.calculate_Mcrd_xx(input)