import numpy as np
import Geometry as Geometry
import Arc_length
import Load_control
import DR_implicite as DR
import Post_pocess
import Structures_tests


# DEFINITIONS
#Elem_type = 1 : beam alu
E1 = 70.390e9 # [N/m^2]
d1 = 0.06 # [m]
t1 = 0.002 # [m]
A1 = np.pi*((d1/2)**2 - (d1/2 - t1)**2) # [m^2]
I1 = (np.pi/4)*((d1/2)**4 - (d1/2 - t1)**4) # [m^4]
W_el1 = I1/(d1/2) # [m^3]
v = 0.3 #[/]
G1 = E1/(2*(1+v)) # [N/m^2]  
J1 = 2*I1  # Approximation pour une section circulaire creuse
rho1 = 2700  # Densité en kg/m^3²
fy1 = 160 #[Mpa]

#Elem_type = 2 : cable
E2 = 72.000e9 # [N/m^2]
E2 = 200.000e9 # [N/m^2]
d2 = 0.005 # [m]
A2 = np.pi*(d2/2)**2 # [m^2]
I2 = 0 # [m^4]
W_el2 = 0 # [m^3]
v = 0.3 #[/]
G2 = E2/(2*(1+v)) # [N/m^2]
J2 = 2*I2  # Pas de rigidité en torsion pour un câble
rho2 = 7850  # Densité en kg/m^3
fy2 = 120 #[Mpa]

#Elem_type = 3 : beam bamboo
E3 = 6.000e9 # [N/m^2]
d3 =  0.145# [m]
t3 =  0.015# [m]
A3 = np.pi*((d3/2)**2 - (d3/2 - t3)**2) # [m^2]µ
I3 = (np.pi/4)*((d3/2)**4 - (d3/2 - t3)**4) # [m^4]
W_el3 = I3/(d3/2) # [m^3]
v = 0.3 #[/]
G3 = E3/(2*(1+v)) # [N/m^2]  
J3 = 2*I3 # Approximation pour une section circulaire creuse
rho3 = 625  # Densité en kg/m^3
fy3 = 50 #[Mpa] !!! bois 

#Elem_type = 4 : cable
E4 = 58.000e9 # [N/m^2]
d4 =  0.01# [m]
A4 = np.pi*(d4/2)**2 # [m^2]
I4 = 0 # [m^4]
W_el4 = 0 # [m^3]
v = 0.3 #[/]
G4 = E4/(2*(1+v)) # [N/m^2]
J4 = 2*I4 # Pas de rigidité en torsion pour un câble
rho4 = 0.24/A4  # Densité en kg/m^3
fy4 = 700 #[Mpa]

d5 = 0.02
A5 = np.pi*(d5/2)**2
I5 = np.pi*(d5/2)**4/4
I5 = 0
E5 = 210e9
G5 = E5/(2*(1+v))
J5 = 2*I5
rho5 = 7850
fy5 = 250 #[Mpa]
W_el5 = I5/(d5/2)

Elem_properties = [[A1, E1, I1, G1, J1,rho1,fy1,W_el1], [A2, E2, I2, G2, J2,rho2,fy2,W_el2],[A3, E3, I3, G3, J3,rho3,fy3,W_el3],[A4, E4, I4, G4, J4,rho4,fy4,W_el4],[A5, E5, I5, G5, J5,rho5,fy5,W_el5]] # Propriétés des éléments : [A, E, I, G, J, rho, fy, W_el]
_lambda = 0

n_mesh = 10# Nombre d'éléments par connexion1
PP = False # Prise en compte du poids propre ou non
linear = False # Prise en compte de la non linéarité ou non

Coordonates_base, Connections_base, Constraints, Forces, Pretension = Structures_tests.treillis(3,3,3)
Coord, Connect, Fixed_DoF, Forces, Elem_Types, Pretension, initial_point_index, rotule_start, rotule_end, U = Geometry.create_mesh(n_mesh, Coordonates_base, Connections_base, Constraints, Forces, Pretension, Elem_properties,PP)

Solveur_type = 1 #choix du solveur : 1 pour load control, 2 pour arc length, 3 pour DR implicite

if Solveur_type == 1:
    U, p_bsc = Load_control.solve(Coord, Connect, Fixed_DoF, Forces, Elem_Types, Elem_properties, Pretension, linear, rotule_start, rotule_end, nb_iteration =200 , nb_step=1 )
elif Solveur_type == 2:
    U, p_bsc, _lambda = Arc_length.solve(Coord, Connect, Fixed_DoF, Forces, Elem_Types, Elem_properties, Pretension, linear, rotule_start, rotule_end,U,nb_step = 100)
elif Solveur_type == 3:
    U, p_bsc = DR.solve(Coord, Connect, Fixed_DoF, Forces,Elem_Types, Elem_properties, Pretension, rotule_start, rotule_end, U, max_iter = 500)


#Output
Post_pocess.post_process(Coord, U, Connect, Elem_Types, initial_point_index, _lambda, p_bsc)

