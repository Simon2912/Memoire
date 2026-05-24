import numpy as np
import matplotlib.pyplot as plt

def solve(Coord, Connect, Fixed_DoF, P_ref,Elem_Types, Elem_properties, dlo, linear, rotule_start, rotule_end,U, nb_step):
    No_Ddl = len(Coord[1])*6 
    No_Elem = len(Connect)  

    AE_Elem = np.zeros(No_Elem)
    EI_Elem = np.zeros(No_Elem)
    GJ_Elem = np.zeros(No_Elem)

    #extraction des propriétés des éléments
    for i in range(No_Elem):
        AE_Elem[i] = Elem_properties[int(Elem_Types[i])-1][0] * Elem_properties[int(Elem_Types[i])-1][1]
        EI_Elem[i] = Elem_properties[int(Elem_Types[i])-1][1] * Elem_properties[int(Elem_Types[i])-1][2]
        GJ_Elem[i] = Elem_properties[int(Elem_Types[i])-1][3] * Elem_properties[int(Elem_Types[i])-1][4]
   
    L_Elem = np.zeros(No_Elem)
    r_C = np.zeros((No_Elem,12,12)) 
    Assemblage = np.zeros((No_Elem, 12))

    def cosinus(i):        
        #calcul de l'orientation de l'élément i
        L_x = Coord[0][int(Connect[i][1])] + U[int(Connect[i][1])*6] - Coord[0][int(Connect[i][0])] - U[int(Connect[i][0])*6]
        L_y = Coord[1][int(Connect[i][1])] + U[int(Connect[i][1])*6+1] - Coord[1][int(Connect[i][0])] - U[int(Connect[i][0])*6+1]
        L_z = Coord[2][int(Connect[i][1])] + U[int(Connect[i][1])*6+2] - Coord[2][int(Connect[i][0])] - U[int(Connect[i][0])*6+2]
        L = np.sqrt(L_x**2 + L_y**2 + L_z**2)
       
        if L_Elem[i] == 0:
            print("Erreur : élément de longueur nulle détecté entre les nœuds " + str(Connect[i][0]) + " et " + str(Connect[i][1]) + ". Veuillez vérifier les coordonnées des nœuds.")
            raise ValueError("Longueur d'élément nulle.")
        cx = L_x / L  
        cy = L_y / L  
        cz = L_z / L  
        return cx, cy, cz

    def pretension():
        # Calcul des forces équivalentes dues à la prétension
        equivalent_loads = np.zeros(No_Ddl)
        equivalent_tension = np.zeros(No_Elem)
        for i in range(No_Elem):
            if dlo[i] != 0:
                equivalent_tension[i] = -dlo[i]*AE_Elem[i]/L_Elem[i]
                cx, cy, cz = cosinus(i)
                equivalent_loads[int(Connect[i][0]*6)] += -equivalent_tension[i] * -cx
                equivalent_loads[int(Connect[i][0]*6+1)] += -equivalent_tension[i] * -cy
                equivalent_loads[int(Connect[i][0]*6+2)] += -equivalent_tension[i] * -cz
                equivalent_loads[int(Connect[i][1]*6)] += -equivalent_tension[i] * cx
                equivalent_loads[int(Connect[i][1]*6+1)] += -equivalent_tension[i] * cy
                equivalent_loads[int(Connect[i][1]*6+2)] += -equivalent_tension[i] * cz
        return equivalent_tension, equivalent_loads

    for i in range(No_Elem) : 
        x1 = Coord[0][int(Connect[i][0])]
        y1 = Coord[1][int(Connect[i][0])]
        z1 = Coord[2][int(Connect[i][0])]
        x2 = Coord[0][int(Connect[i][1])]
        y2 = Coord[1][int(Connect[i][1])]
        z2 = Coord[2][int(Connect[i][1])]
        #evaluation de la matrice de rotation r_C (local vers global) pour chaque élément
        n1 = np.array([x1, y1, z1])
        n2 = np.array([x2, y2, z2])
        v_ref=np.array([0,0,1])
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
      
        if L_Elem[i] == 0:
            print("Erreur : élément de longueur nulle détecté entre les nœuds " + str(Connect[i][0]) + " et " + str(Connect[i][1]) + ". Veuillez vérifier les coordonnées des nœuds.")
            raise ValueError("Longueur d'élément nulle.")
       
       
        # Auxiliary matrices for the assembly: positioning of local matrices in the global matrix
        Assemblage[i] = np.array([int(Connect[i][0]*6), int(Connect[i][0]*6+1), int(Connect[i][0]*6+2), int(Connect[i][0]*6+3), int(Connect[i][0]*6+4), int(Connect[i][0]*6+5),
                                int(Connect[i][1]*6), int(Connect[i][1]*6+1), int(Connect[i][1]*6+2), int(Connect[i][1]*6+3), int(Connect[i][1]*6+4), int(Connect[i][1]*6+5)])
        Assemblage = Assemblage.astype(int)

    #initialisation des vecteurs et matrices
    u_loc = np.zeros((No_Elem, 12))  
    u_glob = np.zeros((No_Elem, 12))  
    u_bsc = np.zeros((No_Elem, 6))  
    p_bsc = np.zeros((No_Elem, 6))  
    Compatibility_matrix_local_basic_non_linear = np.zeros((No_Elem, 6, 12)) 
    G1 = np.zeros((No_Elem,12,12))  
    k_mat = np.zeros((No_Elem,12,12))
    k_geom = np.zeros((No_Elem,12,12))
    k_loc = np.zeros((No_Elem,12,12))
    k_glob = np.zeros((No_Elem,12,12))
    p_loc = np.zeros((No_Elem, 12)) 
    p_global = np.zeros((No_Elem, 12))
    U = np.zeros(No_Ddl)
    F = np.zeros(No_Ddl)
    _lambda = 0.0
    step = 0
    l0 = 1/nb_step

    while _lambda < 1.0 and step < 10*nb_step:
        step += 1
        equivalent_tension, equivalent_loads = pretension()
       
        P = P_ref + equivalent_loads
         
        for i in range(No_Elem) :
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
            l = np.sqrt((L_Elem + delta_ulx) ** 2 + delta_uly ** 2 + delta_ulz ** 2)

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
            
            # Matrice de compatibilité (local vers basic non linéaire)
            A_i = np.array([
                [-cB*cG,  -sB,   -sG,   0,       0,   0,     cB*cG,  sB,     sG,     0,      0,    0],
                [-sB/li,   cB/li,  0,     sG,      0,   cG,    sB/li,   -cB/li,  0,      0,      0,    0],
                [-sB/li,   cB/li,  0,     0,       0,   0,     sB/li,   -cB/li,  0,      sG,     0,    cG],
                [sG/li,   0,     -cG/li,  -sB,     cB,  0,     -sG/li,   0,      cG/li,  0,      0,    0],
                [sG/li,   0,     -cG/li,  0,       0,   0,     -sG/li,   0,      cG/li,  -sB,    cB,   0],
                [0,       0,     0,      -cB*cG,  -sB,  sG,    0,     0,      0,      cB*cG,  sB,   -sG]
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
                u_bsc[i][3] = u_loc[i][4] + gamma[i]
                u_bsc[i][4] = u_loc[i][10] + gamma[i]
                u_bsc[i][5] = u_loc[i][9]- u_loc[i][3]

            # Application des conditions de rotule
            if rotule_start[i] == 1:
                k_bsc[1, :, i] = 0; k_bsc[:, 1, i] = 0; k_bsc[1, 1, i] = 1e-8  
                k_bsc[3, :, i] = 0; k_bsc[:, 3, i] = 0; k_bsc[3, 3, i] = 1e-8 
            if rotule_end[i] == 1:
                k_bsc[2, :, i] = 0; k_bsc[:, 2, i] = 0; k_bsc[2, 2, i] = 1e-8
                k_bsc[4, :, i] = 0; k_bsc[:, 4, i] = 0; k_bsc[4, 4, i] = 1e-8

            p_bsc[i] = k_bsc[:,:,i] @ u_bsc[i] + np.array([equivalent_tension[i], 0, 0,0,0,0])
            
            #Vérification de la non-compression des câbles
            if p_bsc[i][0] + equivalent_tension[i] < 0: 
                    if Elem_Types[i] == 2 or Elem_Types[i] == 4:
                        p_bsc[i][0] = 0
                        AE_Elem[i] = 1e-8  
                        k_bsc[0,0,i] = AE_Elem[i] / L_Elem[i]
            elif Elem_Types[i] == 2 or Elem_Types[i] == 4:  
                AE_Elem[i] = Elem_properties[int(Elem_Types[i])-1][0] * Elem_properties[int(Elem_Types[i])-1][1]
                k_bsc[0,0,i] = AE_Elem[i] / L_Elem[i]

            p_loc[i] = np.transpose(Compatibility_matrix_local_basic_non_linear[i]) @ p_bsc[i]
            p_global[i] = np.transpose(r_C[i]) @ p_loc[i]
            P_r[Assemblage[i]] += p_global[i]
            
        for i in range(No_Elem) : 
            cB = cosB[i]
            sB = sinB[i]
            cG = cosG[i]
            sG = sinG[i]
            L = l[i]
            
        
            Gi = 1/L*np.array([
    [ sB**2 + sG**2,  -sB*cB,   sG*cG,      0, 0, 0,   -sB**2 - sG**2,  sB*cB + 0,   0 + sG*cG,   0, 0, 0],
    [-sB*cB + 0,      cB**2,    0,          0, 0, 0,   cB*sB + 0,       -cB**2 + 0,  0 + 0,       0, 0, 0],
    [ 0 - sG*cG,      0,        0 + cG**2,  0, 0, 0,   0 + sG*cG,       0 + 0,       0 - cG**2,   0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-sB**2 - sG**2,  sB*cB,     -sG*cG,    0, 0, 0,   sB**2 + sG**2,  -sB*cB + 0,   0 - sG*cG,   0, 0, 0],
    [ sB*cB,          -cB**2,    0,         0, 0, 0,   -cB*sB + 0,     cB**2 + 0,    0 + 0,       0, 0, 0],
    [ 0 + sG*cG,      0,         -cG**2,    0, 0, 0,   0 - sG*cG,      0 + 0,        0 + cG**2,   0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ])
            G1[i] = Gi

        K_str = np.zeros((No_Ddl,No_Ddl))

        for i in range(No_Elem) : 
            k_mat[i] = np.transpose(Compatibility_matrix_local_basic_non_linear[i]) @ k_bsc[:,:,i] @ Compatibility_matrix_local_basic_non_linear[i]
            if linear == True:
                k_loc[i] = k_mat[i]

            else:
                k_geom[i] = p_bsc[i,0]*G1[i] 
                k_loc[i] = k_mat[i] + k_geom[i]
            k_glob[i] = np.transpose(r_C[i]) @ k_loc[i] @ r_C[i]
            for j in range(len(k_glob[i])) : 
                for k in range(len(k_glob[i])) :
                    K_str[Assemblage[i][j]][Assemblage[i][k]] += k_glob[i][j][k]
       
        K_str += np.eye(K_str.shape[0]) * 1e-7 # ajout d'une petite rigidité pour éviter les problèmes de singularité

        #application des conditions de Dirichlet
        for i in range(len(Fixed_DoF)):
            K_str[Fixed_DoF[i],:] = 0
            K_str[:,Fixed_DoF[i]] = 0
            K_str[Fixed_DoF[i],Fixed_DoF[i]] = 1
            P[Fixed_DoF[i]] = 0

        delta_U = np.linalg.solve(K_str,P)
        
        for i in range(No_Elem) :
            u_glob[i] = delta_U[Assemblage[i]]
            u_loc[i] = r_C[i] @ u_glob[i]
            
        norm_square = np.dot(delta_U, delta_U)
        d_lambda = l0/np.sqrt(1+norm_square) * np.sign(np.dot(P.transpose(), delta_U))
        if _lambda + d_lambda > 1.0:
            d_lambda = 1.0 - _lambda
        U += d_lambda * delta_U
        F += d_lambda * P
        N = p_bsc[0,0]
        _lambda += d_lambda
        

    print("lambda final :", _lambda,"at step:",step)
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
        
    for i in range(No_Elem) :
        
        u_bsc[i][0] = l[i] - L_Elem[i]
        u_bsc[i][1] = u_loc[i][5] - beta[i]
        u_bsc[i][2] = u_loc[i][11] - beta[i]
        u_bsc[i][3] = u_loc[i][4] + gamma[i]
        u_bsc[i][4] = u_loc[i][10] + gamma[i]
        u_bsc[i][5] = u_loc[i][9]- u_loc[i][3]
        
        p_bsc[i] = k_bsc[:,:,i] @ u_bsc[i]
    
    for i in range(No_Elem):
        p_bsc[i] += np.array([equivalent_tension[i], 0, 0,0,0,0])
    return U, p_bsc, _lambda
           