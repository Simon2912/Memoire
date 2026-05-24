import numpy as np
import Post_process as Post_process
import Load_control
import Structures_tests
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import Structures_tests
import DR_truss as DR
import Geometry 
# import Post_process # Commenté pour que le code soit autonome ici
def K_reduced(N,e0,L0,beta,Elem_properties,e,PP,convergence,Elem_type):
    n_mesh = 10
    M = N *e
    Coordonates_base, Connections_base, Constraints, Forces, Pretension = Structures_tests.poutres(e0,L0,int(N),int(M),Elem_type)
    Coord, Connect, Fixed_DoF, Forces, Elem_Types, Pretension, initial_point_index = Geometry.create_mesh_beam(n_mesh, Coordonates_base, Connections_base, Constraints, Forces, Pretension,Elem_properties,PP,beta)
    U, p_bsc,convergence = Load_control.solve(Coord, Connect, Fixed_DoF, Forces, Elem_Types, Elem_properties ,convergence, linear =False, nb_iteration=200)
    
    N_ed = max(p_bsc[:,0])
    M_ed1 = max(max(abs(p_bsc[:,1])),max(abs(p_bsc[:,2])))
    M_ed2 = max(max(abs(p_bsc[:,3])),max(abs(p_bsc[:,4])))
    E = Elem_properties[Elem_type-1][1]
    A = Elem_properties[Elem_type-1][0]
    dx_nonlineaire = U[-6]
    dx_linéair = N/(E*A/L0)
    K_ratio = dx_linéair/dx_nonlineaire
    if K_ratio >1:
        K_ratio = 1
    return K_ratio, convergence, N_ed, M_ed1, M_ed2


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
rho1 = 2700  # Densité en kg/m^3
fy1 = 160 #[Mpa]

#Elem_type = 2 : cable
E2 = 72.000e9 # [N/m^2]
d2 = 0.008 # [m]
A2 = np.pi*(d2/2)**2 # [m^2]
I2 = 0 # [m^4]
W_el2 = 0 # [m^3]
v = 0.3 #[/]
G2 = E2/(2*(1+v)) # [N/m^2]
J2 = 2*I2  # Pas de rigidité en torsion pour un câble
rho2 = 0  # Densité en kg/m^3
fy2 = 120 #[Mpa]

#Elem_type = 3 : beam bamboo
E3 = 6.000e9 # [N/m^2]
d3 =  0.145# [m]
t3 =  0.015# [m]
A3 = np.pi*((d3/2)**2 - (d3/2 - t3)**2) # [m^2]
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
Elem_properties = [[A1, E1, I1, G1, J1,rho1,fy1,W_el1], [A2, E2, I2, G2, J2,rho2,fy2,W_el2],[A3, E3, I3, G3, J3,rho3,fy3,W_el3],[A4, E4, I4, G4, J4,rho4,fy4,W_el4]]
F = 300 #kg
P = F * 9.81 #N
Coordonates_base, Connections, Constraints, Forces, Pretension = Structures_tests.simplex(2000)
Fixed_Dof_Indices, Forces_Vec, Elem_types, initial_point_index = Geometry.create_mesh_truss(Coordonates_base,Connections,Constraints, Forces)

# CORRECTION CRITIQUE : Transposition des coordonnées (.T)
# Le solveur veut (3, N), trois_cables donne (N, 3)
Coord = Coordonates_base.T 

non_linéarité = True# True pour activer la boucle de non-linéarité, False pour une seule itération linéaire
e0 =300 #Imperfection initiale (0 pour pas d'imperfection, sinon une valeur comme 100 pour e0/L0 = 1/100 d'imperfection)
e = 0.00  #Exentricité de la charge en mètre
PP = 0 #prise en compte du poids propre ou non
Efforts_internes = np.zeros((len(Connections),3))
L0 = np.zeros(len(Connections))
for i in range(len(Connections)):
    n1, n2 = int(Connections[i][0]), int(Connections[i][1])
    vec0 = Coordonates_base[n2] - Coordonates_base[n1]
    L0[i] = np.linalg.norm(vec0)

current_axial_forces = np.zeros(len(Connections))
K_ratio = np.ones(len(Connections)) 
iter = 0
iter_max = 20

convergence = True
while non_linéarité: 
    iter +=1
    U, axial_forces,convergence = DR.solve_truss_3d_DR(Coord, Connections, Fixed_Dof_Indices, Forces_Vec, Elem_types, Elem_properties, Pretension,K_ratio,convergence)
    if np.linalg.norm(axial_forces - current_axial_forces) < 1:
        print("Convergence des forces axiales atteinte.") 
        break

    current_axial_forces = axial_forces.copy()
    for i in range(len(axial_forces)):
        if Elem_types[i] == 2 or Elem_types[i] == 4:
            Efforts_internes[i][0] = axial_forces[i]
    for i in range(len(Connections)):
        if Elem_types[i] == 1 or Elem_types[i] == 3:
            
            n1, n2 = int(Connections[i][0]), int(Connections[i][1])
            disp_n1 = U[3*n1 : 3*n1+3] 
            disp_n2 = U[3*n2 : 3*n2+3]
            vec0 = (Coordonates_base[n2] + disp_n2) - (Coordonates_base[n1] + disp_n1)
            longueur = np.linalg.norm(vec0)
            dz = vec0[2] 
            beta = abs(np.degrees(np.arcsin(dz / longueur))) # angle entre la barre et le sol pour appliquer le poids propre
            K_ratio_previous = K_ratio[i]
            K_ratio[i],convergence, N_ed, M_ed1, M_ed2 = K_reduced(axial_forces[i],e0,L0[i],beta,Elem_properties,e,PP,convergence,Elem_types[i])
            K_ratio[i] = (K_ratio[i] + K_ratio_previous)/2
            Efforts_internes[i][0] = N_ed
            Efforts_internes[i][1] = M_ed1
            Efforts_internes[i][2] = M_ed2

        else:
            K_ratio[i] = 1.0
    if iter >=iter_max:
        print("Nombre maximum d'itérations atteint sans convergence.")
        break

if non_linéarité == False:
    U, axial_forces,convergence = DR.solve_truss_3d_DR(Coord, Connections, Fixed_Dof_Indices, Forces_Vec, Elem_types, Elem_properties, Pretension,K_ratio,convergence)
for i in range(len(axial_forces)):
    print(f"Element {i} : force axiale = {axial_forces[i]:.3f} N, K_ratio = {K_ratio[i]:.3f}")
Post_process.post_process(Coord,U, Connections[:,:2], Elem_types, initial_point_index, convergence)