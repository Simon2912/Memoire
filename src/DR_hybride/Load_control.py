import numpy as np
import matplotlib.pyplot as plt

def solve(Coord, Connect, Fixed_DoF, P,Elem_Types, Elem_properties, convergence,linear, nb_iteration,):
    No_Ddl = len(Coord[1])*6  # 3 DoF per node
    No_Elem = len(Connect)  # Number of elements
    AE_Elem = np.zeros(No_Elem)
    EI_Elem = np.zeros(No_Elem)
    GJ_Elem = np.zeros(No_Elem)
   
    for i in range(No_Elem):
        AE_Elem[i] = Elem_properties[int(Elem_Types[i])-1][0] * Elem_properties[int(Elem_Types[i])-1][1]
        EI_Elem[i] = Elem_properties[int(Elem_Types[i])-1][1] * Elem_properties[int(Elem_Types[i])-1][2]
        GJ_Elem[i] = Elem_properties[int(Elem_Types[i])-1][3] * Elem_properties[int(Elem_Types[i])-1][4]
    
    L_Elem = np.zeros(No_Elem)
    r_C = np.zeros((No_Elem,12,12)) #MODIFIE
    Assemblage = np.zeros((No_Elem, 12))

    for i in range(No_Elem) : 
        U = np.zeros(No_Ddl)
        x1 = Coord[0][int(Connect[i][0])]
        y1 = Coord[1][int(Connect[i][0])]
        z1 = Coord[2][int(Connect[i][0])]
        x2 = Coord[0][int(Connect[i][1])]
        y2 = Coord[1][int(Connect[i][1])]
        z2 = Coord[2][int(Connect[i][1])]
        # 1. Element's length
        n1 = np.array([x1, y1, z1])
        n2 = np.array([x2, y2, z2])
        v_ref=np.array([0,0,1])
        # direction locale x
        ex = n2 - n1
        ex = ex / np.linalg.norm(ex)
        ey = np.cross(ex, v_ref)
        
        if np.linalg.norm(ey) < 1e-8:
            v_ref = np.array([0,1,0])
            ey = np.cross(ex, v_ref)
            
        ey = ey / np.linalg.norm(ey)
        ez = np.cross(ey, ex)
        ey*=-1

        t = np.vstack((ex, ey, ez))
        T = np.zeros((12, 12))
        for j in range(4):
            T[j*3:(j+1)*3, j*3:(j+1)*3] = t

        r_C[i] = T

        L_Elem[i] = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        
        k = np.array([[AE_Elem[0]/L_Elem[0], 0, 0, 0, 0, 0,-AE_Elem[0]/L_Elem[0],0, 0, 0, 0, 0],
                  [0, 12*EI_Elem[0]/(L_Elem[0]**3), 0, 0, 0, 6*EI_Elem[0]/(L_Elem[0]**2), 0,-12*EI_Elem[0]/(L_Elem[0]**3), 0, 0, 0, 6*EI_Elem[0]/(L_Elem[0]**2)],
                  [0, 0, 12*EI_Elem[0]/(L_Elem[0]**3), 0, -6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, -12*EI_Elem[0]/(L_Elem[0]**3), 0, -6*EI_Elem[0]/(L_Elem[0]**2),0],
                  [0, 0, 0, GJ_Elem[0]/L_Elem[0], 0, 0, 0, 0, 0, -GJ_Elem[0]/L_Elem[0], 0, 0],
                    [0, 0, -6*EI_Elem[0]/(L_Elem[0]**2), 0, 4*EI_Elem[0]/L_Elem[0], 0, 0, 0, 6*EI_Elem[0]/(L_Elem[0]**2), 0, 2*EI_Elem[0]/L_Elem[0], 0],
                    [0, 6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, 4*EI_Elem[0]/L_Elem[0], 0, -6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, 2*EI_Elem[0]/L_Elem[0]],
                    [-AE_Elem[0]/L_Elem[0], 0, 0, 0, 0, 0, AE_Elem[0]/L_Elem[0],0, 0, 0, 0, 0],
                    [0, -12*EI_Elem[0]/(L_Elem[0]**3), 0, 0, 0, -6*EI_Elem[0]/(L_Elem[0]**2), 0, 12*EI_Elem[0]/(L_Elem[0]**3), 0, 0, 0, -6*EI_Elem[0]/(L_Elem[0]**2)],
                    [0, 0, -12*EI_Elem[0]/(L_Elem[0]**3), 0, 6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, 12*EI_Elem[0]/(L_Elem[0]**3), 0,6*EI_Elem[0]/(L_Elem[0]**2), 0],
                    [0, 0, 0, -GJ_Elem[0]/L_Elem[0], 0, 0, 0, 0, 0, GJ_Elem[0]/L_Elem[0], 0, 0],
                    [0, 0,-6*EI_Elem[0]/(L_Elem[0]**2), 0, 2*EI_Elem[0]/L_Elem[0], 0, 0, 0, 6*EI_Elem[0]/(L_Elem[0]**2), 0, 4*EI_Elem[0]/L_Elem[0], 0],
                    [0, 6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, 2*EI_Elem[0]/L_Elem[0], 0, -6*EI_Elem[0]/(L_Elem[0]**2), 0, 0, 0, 4*EI_Elem[0]/L_Elem[0]]])
        
 
        if L_Elem[i] == 0:
            print("Erreur : élément de longueur nulle détecté entre les nœuds " + str(Connect[i][0]) + " et " + str(Connect[i][1]) + ". Veuillez vérifier les coordonnées des nœuds.")
            raise ValueError("Longueur d'élément nulle.")
              
        Assemblage[i] = np.array([Connect[i][0]*6, 
                                Connect[i][0]*6+1,
                                Connect[i][0]*6+2,
                                Connect[i][0]*6+3,
                                Connect[i][0]*6+4,
                                Connect[i][0]*6+5,
                                Connect[i][1]*6,
                                Connect[i][1]*6+1,
                                Connect[i][1]*6+2,
                                Connect[i][1]*6+3,
                                Connect[i][1]*6+4,
                                Connect[i][1]*6+5])
        Assemblage = Assemblage.astype(int)
    # Initialization of vecteurs
    u_loc = np.zeros((No_Elem, 12))  #MODIFIE
    p_loc = np.zeros((No_Elem, 12))  #MODIFIE
    u_glob = np.zeros((No_Elem, 12))  #MODIFIE
    u_bsc = np.zeros((No_Elem, 6))  #MODIFIE
    p_bsc = np.zeros((No_Elem, 6))  #MODIFIE
    p_global = np.zeros((No_Elem, 12))  #MODIFIE
    Compatibility_matrix_local_basic_non_linear = np.zeros((No_Elem, 6, 12))  #MODIFIE
    G1 = np.zeros((No_Elem,12,12))  #MODIFIE
    k_mat = np.zeros((No_Elem,12,12))
    k_geom = np.zeros((No_Elem,12,12))
    k_loc = np.zeros((No_Elem,12,12))
    k_glob = np.zeros((No_Elem,12,12))

    for iter in range(nb_iteration):
        for i in range(No_Elem): 
            u_glob[i] = U[Assemblage[i]]
            u_loc[i] = r_C[i] @ u_glob[i]
        
        delta_ulx = u_loc[:, 6] - u_loc[:, 0]
        delta_uly = u_loc[:, 7] - u_loc[:, 1]
        delta_ulz = u_loc[:, 8] - u_loc[:, 2]
        
        if linear == True:
            beta = np.zeros(No_Elem)
            gamma = np.zeros(No_Elem)
            l = L_Elem
        else:
            beta = np.arctan2(delta_uly, (L_Elem + delta_ulx))
            gamma = np.arctan2(delta_ulz, (L_Elem + delta_ulx))
            l = np.sqrt((L_Elem + delta_ulx) ** 2 + delta_uly ** 2+ delta_ulz ** 2)
            # Calcul de r_cbsc et de u_bsc
        cosB = np.cos(beta)
        sinB = np.sin(beta)
        cosG = np.cos(gamma)
        sinG = np.sin(gamma)
        for i in range(No_Elem):
            cB = cosB[i]
            sB = sinB[i]
            cG = cosG[i]
            sG = sinG[i]
            li  = l[i]

            A_i = np.array([
                [-cB*cG,  -sB,   -sG,   0,       0,   0,     cB*cG,  sB,     sG,     0,      0,    0],
                [-sB/li,   cB/li,  0,     sG,      0,   cG,    sB/li,   -cB/li,  0,      0,      0,    0],
                [-sB/li,   cB/li,  0,     0,       0,   0,     sB/li,   -cB/li,  0,      sG,     0,    cG],
                [-sG/li,   0,     cG/li,  -sB,     cB,  0,     sG/li,   0,      -cG/li,  0,      0,    0],
                [-sG/li,   0,     cG/li,  0,       0,   0,     sG/li,   0,      -cG/li,  -sB,    cB,   0],
                [0,       0,     0,      cB*cG,  sB,  -sG,    0,     0,      0,      cB*cG,  sB,   -sG]
            ])

            Compatibility_matrix_local_basic_non_linear[i] = A_i
        for i in range(No_Elem) : 
            u_bsc[i] = Compatibility_matrix_local_basic_non_linear[i] @ u_loc[i]

        k_bsc = np.array([[AE_Elem / L_Elem,                  np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem)],
                [np.zeros(No_Elem),         4 * EI_Elem / L_Elem,    2 * EI_Elem / L_Elem,         np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem)],
                [np.zeros(No_Elem),         2 * EI_Elem / L_Elem,    4 * EI_Elem / L_Elem,         np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem)],
                [np.zeros(No_Elem),                 np.zeros(No_Elem),         np.zeros(No_Elem),  4 * EI_Elem / L_Elem,    2 * EI_Elem / L_Elem,         np.zeros(No_Elem)],
                [np.zeros(No_Elem),                 np.zeros(No_Elem),         np.zeros(No_Elem),  2 * EI_Elem / L_Elem,    4 * EI_Elem / L_Elem,         np.zeros(No_Elem)],
                [np.zeros(No_Elem),                 np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem),         np.zeros(No_Elem),   GJ_Elem / L_Elem]
                ])

        # Calcul de p_bsc à P_r
        P_r = np.zeros(No_Ddl)
        for i in range(No_Elem) :
            if linear == False:
                u_bsc[i][0] = l[i] - L_Elem[i]
                u_bsc[i][1] = u_loc[i][5] - beta[i]
                u_bsc[i][2] = u_loc[i][11] - beta[i]
                u_bsc[i][3] = u_loc[i][4] - gamma[i]
                u_bsc[i][4] = u_loc[i][10] - gamma[i]
                u_bsc[i][5] = u_loc[i][9]- u_loc[i][3]
           
            p_bsc[i] = k_bsc[:,:,i] @ u_bsc[i]# + np.array([equivalent_tension[i], 0, 0])
            p_loc[i] = np.transpose(Compatibility_matrix_local_basic_non_linear[i]) @ p_bsc[i]
            p_global[i] = np.transpose(r_C[i]) @ p_loc[i]
            P_r[Assemblage[i]] += p_global[i]
   
        R = P - P_r
        for i in range(No_Elem) : 
            cB = cosB[i]
            sB = sinB[i]
            cG = cosG[i]
            sG = sinG[i]
            L = l[i]
                        
            Gi = 1/L*np.array([
[ sB**2 + sG**2, -sB*cB,sB*cG,       0, 0, 0, -sB**2 - sG**2, sB*cB + 0, 0 + sB*cG, 0, 0, 0],
[-sB*cB + 0,      cB**2,0,           0, 0, 0,  cB*sB + 0, -cB**2 + 0, 0 + 0, 0, 0, 0],
[ 0 - sG*cG,      0,          0 + cG**2,       0, 0, 0, 0 + sG*cG, 0 + 0, 0 - cG**2, 0, 0, 0],
[ 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0],
[-sB**2 - sG**2,cB*sB + 0,sG*cG,       0, 0, 0, sB**2 + sG**2, -sB*cB + 0, 0 - sG*cG, 0, 0, 0],
[ sB*cB,-cB**2,    0,           0, 0, 0, -cB*sB + 0, cB**2 + 0, 0 + 0, 0, 0, 0],
[ 0 + sG*cG,0,cG**2,        0, 0, 0, 0 - sG*cG, 0 + 0, 0 + cG**2, 0, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])
            G1[i] = Gi


        G23 = (1/(l**2))*np.array([[-2*cosB*sinB, cosB**2 - sinB**2, np.zeros(No_Elem), 2*cosB*sinB, sinB**2 - cosB**2, np.zeros(No_Elem)],
                                [cosB**2 - sinB**2, 2*cosB*sinB, np.zeros(No_Elem), sinB**2 - cosB**2, -2*cosB*sinB, np.zeros(No_Elem)],
                                [np.zeros(No_Elem), np.zeros(No_Elem), np.zeros(No_Elem), np.zeros(No_Elem),np.zeros(No_Elem), np.zeros(No_Elem)],
                                [2*cosB*sinB, sinB**2 - cosB**2, np.zeros(No_Elem), -2*cosB*sinB, cosB**2 - sinB**2, np.zeros(No_Elem)],
                                [sinB**2 - cosB**2, -2*cosB*sinB, np.zeros(No_Elem), cosB**2 - sinB**2, 2*cosB*sinB, np.zeros(No_Elem)],
                                [np.zeros(No_Elem), np.zeros(No_Elem), np.zeros(No_Elem), np.zeros(No_Elem),np.zeros(No_Elem), np.zeros(No_Elem)]])
        
        K_str = np.zeros((No_Ddl,No_Ddl))
        for i in range(No_Elem) : 
            k_mat[i] = np.transpose(Compatibility_matrix_local_basic_non_linear[i]) @ k_bsc[:,:,i] @ Compatibility_matrix_local_basic_non_linear[i]
            if linear == True:
                k_loc[i] = k_mat[i]
            else:
                k_geom[i] = p_bsc[i,0]*G1[i] #+ (p_bsc[i,1]+p_bsc[i,2])*G23[:,:,i]
                k_loc[i] = k_mat[i] #+ k_geom[i]
                
            k_glob[i] = np.transpose(r_C[i]) @ k_loc[i] @ r_C[i]
            for j in range(len(k_glob[i])) : 
                for k in range(len(k_glob[i])) :
                    K_str[Assemblage[i][j]][Assemblage[i][k]] += k_glob[i][j][k]

        K_str += np.eye(K_str.shape[0]) * 1e-5
        for i in range(len(Fixed_DoF)):
            K_str[Fixed_DoF[i],:] = 0
            K_str[:,Fixed_DoF[i]] = 0
            K_str[Fixed_DoF[i],Fixed_DoF[i]] = 1
            R[Fixed_DoF[i]] = 0
        Residual = np.linalg.norm(R)
                        
        tol_force = 1e-5
        if Residual <= tol_force:
            print("Convergence reached at iteration " + str(iter))
            break
        
        elif iter == nb_iteration - 1:
            print("Convergence not reached with a residual of " + str(Residual))
            convergence = False
            break
        
        U += np.linalg.solve(K_str,R)

    return U, p_bsc,convergence