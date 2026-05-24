import numpy as np
import math
import matplotlib.pyplot as plt

#transfomation d'un vecteur de rotation en quaternion
def rotvec_to_quat(rot_vec):
    angle = np.linalg.norm(rot_vec)
    if angle < 1e-10:
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = rot_vec / angle
    sin_a = math.sin(angle / 2.0)
    return np.array([math.cos(angle / 2.0), axis[0]*sin_a, axis[1]*sin_a, axis[2]*sin_a])

#addition de deux quaternions
def quat_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    q3 = np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])
    return q3

#transformation d'un quaternion en matrice de rotation
def quat_to_matrix(q):
    w, x, y, z = q
    return np.array([
        [1 - 2*y**2 - 2*z**2,     2*x*y - 2*z*w,     2*x*z + 2*y*w],
        [    2*x*y + 2*z*w, 1 - 2*x**2 - 2*z**2,     2*y*z - 2*x*w],
        [    2*x*z - 2*y*w,     2*y*z + 2*x*w, 1 - 2*x**2 - 2*y**2]
    ])

# Extraction du vecteur de rotation d'une matrice de rotation
def extract_spin(R_mat):
    trace = R_mat[0,0] + R_mat[1,1] + R_mat[2,2]
    cos_theta = np.clip((trace - 1.0) / 2.0, -1.0, 1.0)
    angle = math.acos(cos_theta)
    factor = angle / (2.0 * math.sin(angle)) if angle > 1e-6 else 0.5
    rx = factor * (R_mat[2,1] - R_mat[1,2])
    ry = factor * (R_mat[0,2] - R_mat[2,0])
    rz = factor * (R_mat[1,0] - R_mat[0,1])
    return rx, ry, rz

def solve(Coord, Connect, Fixed_DoF, P, Elem_Types, Elem_properties, dlo, rotule_start, rotule_end,U, max_iter):
    omega_star = 1
    a = 1.0
    dt = 1.3/ omega_star
    xi_star = 1/ np.sqrt(2.0 * a)
    amortissement = dt * xi_star * omega_star 
    lambda_v = (1.0 - amortissement) / (1.0 + amortissement)
    print("lambda_v (relaxation sur la vitesse) :", lambda_v)
    lambda_u = (omega_star**2 * dt) / (1.0 + amortissement)
    print("lambda_u (relaxation sur le déplacement) :", lambda_u)
      
    No_Nodes = len(Coord[1])
    No_Ddl = No_Nodes * 6 
    No_Elem = len(Connect)  
    V_prev = np.zeros(No_Ddl)  
    AE_Elem = np.zeros(No_Elem)
    EI_Elem = np.zeros(No_Elem)
    GJ_Elem = np.zeros(No_Elem)

    # Extraction des propriétés des éléments
    for i in range(No_Elem):
        AE_Elem[i] = Elem_properties[int(Elem_Types[i])-1][0] * Elem_properties[int(Elem_Types[i])-1][1]
        EI_Elem[i] = Elem_properties[int(Elem_Types[i])-1][1] * Elem_properties[int(Elem_Types[i])-1][2]
        GJ_Elem[i] = Elem_properties[int(Elem_Types[i])-1][3] * Elem_properties[int(Elem_Types[i])-1][4]
   
    L_Elem = np.zeros(No_Elem)
    r_C = np.zeros((No_Elem,12,12)) 
    Assemblage = np.zeros((No_Elem, 12), dtype=int)
    
    v1_0_list = np.zeros((No_Elem, 3))
    T_0_list = np.zeros((No_Elem, 3, 3))
    equivalent_tension = np.zeros(No_Elem)

    for i in range(No_Elem): 
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
        L_Elem[i] = np.linalg.norm(ex)
        if L_Elem[i] == 0:
            print("Erreur : élément de longueur nulle détecté entre les nœuds " + str(Connect[i][0]) + " et " + str(Connect[i][1]) + ". Veuillez vérifier les coordonnées des nœuds.")
            raise ValueError("Longueur d'élément nulle.")
              
        ex = ex / L_Elem[i]
        ey = np.cross(ex, v_ref)
        
        if np.linalg.norm(ey) < 1e-8:
            v_ref = np.array([0,1,0])
            ey = np.cross(ex, v_ref)
            
        ey = ey / np.linalg.norm(ey)
        ez = np.cross(ey, ex)
        ey*=-1

        t = np.vstack((ex, ey, ez))
        T_0_list[i] = t
        v1_0_list[i] = ex
        
        if L_Elem[i] == 0:
            print("Erreur : élément de longueur nulle détecté entre les nœuds " + str(Connect[i][0]) + " et " + str(Connect[i][1]) + ". Veuillez vérifier les coordonnées des nœuds.")
            raise ValueError("Longueur d'élément nulle.")
        
        Assemblage[i] = [Connect[i][0]*6, Connect[i][0]*6+1, Connect[i][0]*6+2, Connect[i][0]*6+3, Connect[i][0]*6+4, Connect[i][0]*6+5,
                         Connect[i][1]*6, Connect[i][1]*6+1, Connect[i][1]*6+2, Connect[i][1]*6+3, Connect[i][1]*6+4, Connect[i][1]*6+5]
                         
        # Calcul de la prétension pure
        equivalent_tension[i] = -dlo[i] * (AE_Elem[i] / L_Elem[i])

    # Initialisation des variables globales
    Node_Quats = np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (No_Nodes, 1))
    
    # Allocation des matrices de stockage
    u_bsc = np.zeros((No_Elem, 6))  
    p_bsc = np.zeros((No_Elem, 6))  
    p_loc = np.zeros((No_Elem, 12))  
    p_global = np.zeros((No_Elem, 12))  
    k_mat = np.zeros((No_Elem,12,12))
    k_geom = np.zeros((No_Elem,12,12))
    k_loc = np.zeros((No_Elem,12,12))
    k_glob = np.zeros((No_Elem,12,12))
      
    for iter in range(max_iter):
        P_r = np.zeros(No_Ddl)
        K_str = np.zeros((No_Ddl, No_Ddl))
        
        # Conversion des quaternions des noeuds en matrices
        R_noeuds_globaux = np.zeros((No_Nodes, 3, 3))
        for n in range(No_Nodes):
            R_noeuds_globaux[n] = quat_to_matrix(Node_Quats[n])
            
        for i in range(No_Elem):
            node1, node2 = int(Connect[i][0]), int(Connect[i][1])
            u_el = U[Assemblage[i]]
            
            x1 = Coord[0][node1] + u_el[0]; y1 = Coord[1][node1] + u_el[1]; z1 = Coord[2][node1] + u_el[2]
            x2 = Coord[0][node2] + u_el[6]; y2 = Coord[1][node2] + u_el[7]; z2 = Coord[2][node2] + u_el[8]
            
            dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
            Ln = math.sqrt(dx**2 + dy**2 + dz**2)
            v1 = np.array([dx, dy, dz]) / Ln
            
            v1_0 = v1_0_list[i]
            T_0 = T_0_list[i]
            
            cross_v = np.cross(v1_0, v1)
            dot_v = np.clip(np.dot(v1_0, v1), -1.0, 1.0)
            
            if dot_v < -0.999999: 
                axis = np.cross(v1_0, np.array([1.0, 0.0, 0.0]))
                if np.linalg.norm(axis) < 1e-5: 
                    axis = np.cross(v1_0, np.array([0.0, 1.0, 0.0]))
                axis /= np.linalg.norm(axis)
                ax, ay, az = axis
                R_chord = np.array([
                    [1 - 2*(ay**2 + az**2), 2*ax*ay, 2*ax*az],
                    [2*ax*ay, 1 - 2*(ax**2 + az**2), 2*ay*az],
                    [2*ax*az, 2*ay*az, 1 - 2*(ax**2 + ay**2)]
                ])
            else: 
                skew_v = np.array([[0, -cross_v[2], cross_v[1]], 
                                   [cross_v[2], 0, -cross_v[0]], 
                                   [-cross_v[1], cross_v[0], 0]])
                
                R_chord = np.eye(3) + skew_v + (skew_v @ skew_v) / (1.0 + dot_v)
                
            T_corde = T_0 @ R_chord.T
            
            T = np.zeros((12, 12))
            for k in range(4):
                T[k*3:(k+1)*3, k*3:(k+1)*3] = T_corde
            r_C[i] = T
            
            # Extraction des angles purs
            R_noeud1 = R_noeuds_globaux[node1]
            R_noeud2 = R_noeuds_globaux[node2]
            T_0 = T_0_list[i]
            
            R_loc1 = T_corde @ R_noeud1 @ T_0.T
            R_loc2 = T_corde @ R_noeud2 @ T_0.T
            
            theta_loc_1x, theta_loc_1y, theta_loc_1z = extract_spin(R_loc1)
            theta_loc_2x, theta_loc_2y, theta_loc_2z = extract_spin(R_loc2)
            
            # Formulation Basique
            L = L_Elem[i]
            AE, EI, GJ = AE_Elem[i], EI_Elem[i], GJ_Elem[i]

            u_bsc[i][0] = Ln - L
            u_bsc[i][1] = theta_loc_1y
            u_bsc[i][2] = theta_loc_2y
            u_bsc[i][3] = theta_loc_1z
            u_bsc[i][4] = theta_loc_2z
            u_bsc[i][5] = theta_loc_2x - theta_loc_1x
            
            # Matrice de rigidité basique
            k_bsc_mat = np.zeros((6, 6))
            k_bsc_mat[0, 0] = AE / L
            k_bsc_mat[1, 1] = 4 * EI / L; k_bsc_mat[1, 2] = 2 * EI / L
            k_bsc_mat[2, 1] = 2 * EI / L; k_bsc_mat[2, 2] = 4 * EI / L
            k_bsc_mat[3, 3] = 4 * EI / L; k_bsc_mat[3, 4] = 2 * EI / L
            k_bsc_mat[4, 3] = 2 * EI / L; k_bsc_mat[4, 4] = 4 * EI / L
            k_bsc_mat[5, 5] = GJ / L

            # Application des conditions de rotule
            if rotule_start[i] == 1:
                k_bsc_mat[1, :] = 0; k_bsc_mat[:, 1] = 0; k_bsc_mat[1, 1] = 1e-8
                k_bsc_mat[3, :] = 0; k_bsc_mat[:, 3] = 0; k_bsc_mat[3, 3] = 1e-8
            if rotule_end[i] == 1:
                k_bsc_mat[2, :] = 0; k_bsc_mat[:, 2] = 0; k_bsc_mat[2, 2] = 1e-8
                k_bsc_mat[4, :] = 0; k_bsc_mat[:, 4] = 0; k_bsc_mat[4, 4] = 1e-8
                
            p_bsc[i] = k_bsc_mat @ u_bsc[i]
            
            # application de la prétension
            p_bsc[i][0] += equivalent_tension[i]
            
            # Vérification des câbles
            if int(Elem_Types[i]) in [2, 4]:
                if p_bsc[i][0] < 0:
                    p_bsc[i][0] = 0.0 # Force de compression annulée
            
            B = np.zeros((12, 6))
            B[0, 0] = -1.0; B[6, 0] = 1.0; B[3, 5] = -1.0; B[9, 5] = 1.0; B[4, 1] = 1.0;    B[10, 2] = 1.0; B[2, 1] = -1.0/L
            B[2, 2] = -1.0/L; B[8, 1] = 1.0/L;  B[8, 2] = 1.0/L; B[5, 3] = 1.0;    B[11, 4] = 1.0; B[1, 3] = 1.0/L;  B[1, 4] = 1.0/L; B[7, 3] = -1.0/L; B[7, 4] = -1.0/L
            
            p_loc[i] = B @ p_bsc[i]
            
            k_mat[i] = B @ k_bsc_mat @ B.T
            
            # Matrice géométrique K_G
            F_x2 = p_loc[i][6]; M_x2 = p_loc[i][9]
            M_y1 = p_loc[i][4]; M_z1 = p_loc[i][5]
            M_y2 = p_loc[i][10]; M_z2 = p_loc[i][11]
            
            A_el = Elem_properties[int(Elem_Types[i]) - 1][0]
            I_el = Elem_properties[int(Elem_Types[i]) - 1][2]
            Ip_A = (2.0 * I_el) / A_el 
            
            K_G = np.zeros((12, 12))
            K_G[0, 0] = F_x2 / L; K_G[0, 6] = -F_x2 / L
            K_G[1, 1] = 6.0 * F_x2 / (5.0 * L); K_G[1, 3] = M_y1 / L; K_G[1, 4] = M_x2 / L; K_G[1, 5] = F_x2 / 10.0
            K_G[1, 7] = -6.0 * F_x2 / (5.0 * L); K_G[1, 9] = M_y2 / L; K_G[1, 10] = -M_x2 / L; K_G[1, 11] = F_x2 / 10.0
            K_G[2, 2] = 6.0 * F_x2 / (5.0 * L); K_G[2, 3] = M_z1 / L; K_G[2, 4] = -F_x2 / 10.0; K_G[2, 5] = M_x2 / L
            K_G[2, 8] = -6.0 * F_x2 / (5.0 * L); K_G[2, 9] = M_z2 / L; K_G[2, 10] = -F_x2 / 10.0; K_G[2, 11] = -M_x2 / L
            K_G[3, 3] = (F_x2 * Ip_A) / L; K_G[3, 4] = -(2.0 * M_z1 - M_z2) / 6.0; K_G[3, 5] = (2.0 * M_y1 - M_y2) / 6.0; K_G[3, 7] = -M_y1 / L
            K_G[3, 8] = -M_z1 / L; K_G[3, 9] = -(F_x2 * Ip_A) / L; K_G[3, 10] = -(M_z1 + M_z2) / 6.0; K_G[3, 11] = (M_y1 + M_y2) / 6.0
            K_G[4, 4] = 2.0 * F_x2 * L / 15.0; K_G[4, 7] = -M_x2 / L; K_G[4, 8] = F_x2 / 10.0 
            K_G[4, 9] = -(M_z1 + M_z2) / 6.0; K_G[4, 10] = -F_x2 * L / 30.0; K_G[4, 11] = M_x2 / 2.0
            K_G[5, 5] = 2.0 * F_x2 * L / 15.0; K_G[5, 7] = -F_x2 / 10.0; K_G[5, 8] = -M_x2 / L
            K_G[5, 9] = (M_y1 + M_y2) / 6.0; K_G[5, 10] = -M_x2 / 2.0; K_G[5, 11] = -F_x2 * L / 30.0
            K_G[6, 6] = F_x2 / L
            K_G[7, 7] = 6.0 * F_x2 / (5.0 * L); K_G[7, 9] = -M_y2 / L; K_G[7, 10] = M_x2 / L; K_G[7, 11] = -F_x2 / 10.0
            K_G[8, 8] = 6.0 * F_x2 / (5.0 * L); K_G[8, 9] = -M_z2 / L; K_G[8, 10] = F_x2 / 10.0; K_G[8, 11] = M_x2 / L
            K_G[9, 9] = (F_x2 * Ip_A) / L; K_G[9, 10] = (M_z1 - 2.0 * M_z2) / 6.0; K_G[9, 11] = -(M_y1 - 2.0 * M_y2) / 6.0
            K_G[10, 10] = 2.0 * F_x2 * L / 15.0; K_G[11, 11] = 2.0 * F_x2 * L / 15.0
            
            k_geom[i] = K_G + K_G.T - np.diag(K_G.diagonal())
            
            k_loc[i] = k_mat[i] + k_geom[i]
            k_glob[i] = r_C[i].T @ k_loc[i] @ r_C[i]
            p_global[i] = r_C[i].T @ p_loc[i]
            
            # Assemblage
            for j in range(12):
                P_r[Assemblage[i][j]] += p_global[i][j]
                for k in range(12):
                    K_str[Assemblage[i][j], Assemblage[i][k]] += k_glob[i][j, k]

        K_str += np.eye(K_str.shape[0]) * 1e-6 # ajout d'une petite rigidité pour éviter les problèmes de singularité
        R = P - P_r

        #application des conditions de Dirichlet
        for dof in Fixed_DoF:
            K_str[dof, :] = 0.0
            K_str[:, dof] = 0.0
            K_str[dof, dof] = 1.0
            R[dof] = 0.0
            
        Residual = np.linalg.norm(R)
        
        # Vérification de la convergence
        tol_force = 1e-5
        if Residual <= tol_force:
            print(f"  -> Convergence atteinte à l'itération {iter} (Résidu: {Residual:.2e})")
            for i in range(No_Elem):
                # ajout d'une petite rigidité pour éviter les problèmes de singularité
                if (Elem_Types[i] in [2, 4]) and (p_bsc[i,0] == 0.0): 
                    print(f"     Warning: Câble {i} est détendu (compression annulée).")
            break
        elif iter == max_iter - 1:
            print(f"  -> NON CONVERGENCE : Max itérations atteintes. Résidu: {Residual:.2e}")
            
        du_r = np.linalg.solve(K_str, R) 
        delta_U = dt * (lambda_v * V_prev + lambda_u * du_r)
        U += delta_U

        # Mise à jour des quaternions
        for n in range(No_Nodes):
            idx = n * 6
            delta_rot = delta_U[idx+3:idx+6]
            dq = rotvec_to_quat(delta_rot)
            Node_Quats[n] = quat_mult(dq, Node_Quats[n])

        V_prev = delta_U / dt
    
   
    
    return U, p_bsc