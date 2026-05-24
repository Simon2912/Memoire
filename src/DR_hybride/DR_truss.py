import numpy as np
import matplotlib.pyplot as plt


def solve_truss_3d_DR(Coord, Connect, Fixed_DoF, Forces, Elem_Types, Elem_properties, dlo,K_ratio,convergence):
    # --- 1. Initialisation ---
    # CORRECTION : On s'assure que Coord est bien (3, N)
    # Si l'utilisateur passe (N, 3), Coord.shape[1] donnerait 3 au lieu de N.
    # Ici, on suppose que l'entrée est correcte (3, N) grâce à la correction dans le main.
    
    No_Nodes = Coord.shape[1] 
    No_Ddl = No_Nodes * 3 
    No_Elem = len(Connect)
    
    print(f"--- Structure 3D : {No_Nodes} oeuds, {No_Elem} elements, {No_Ddl} DDL ---")

    # Vecteurs d'état
    U = np.zeros(No_Ddl)        
    V = np.zeros(No_Ddl)        
    KE_list = []
    # Pré-calcul des propriétés élémentaires
    EA_Elem = np.zeros(No_Elem)
    L0_Elem = np.zeros(No_Elem)
    
    for i in range(No_Elem):
        type_idx = int(Elem_Types[i]) - 1
        EA_Elem[i] = Elem_properties[type_idx][0] * Elem_properties[type_idx][1]
        
        n1, n2 = int(Connect[i][0]), int(Connect[i][1])
        # Ici, Coord doit être (3, N) pour que [:, n2] renvoie le vecteur [x,y,z]
        vec0 = Coord[:, n2] - Coord[:, n1]
        L0_Elem[i] = np.linalg.norm(vec0)
    # Paramètres DR
    dt = 0.01                 
    mass_factor = 2       
    tol_residual = 1e-5     
    max_steps = 10000 # Augmenté un peu pour la convergence
    KE_prev = 0.0              
    step = 0
    # --- 2. Boucle Temporelle ---
    while True:
        step += 1
        
        P_int = np.zeros(No_Ddl)
        K_diag_global = np.zeros(No_Ddl) 
        axial_forces = np.zeros(No_Elem)
        equivalent_loads = np.zeros(No_Ddl)
        for i in range(No_Elem):
            n1, n2 = int(Connect[i][0]), int(Connect[i][1])
            idx1 = [n1*3, n1*3+1, n1*3+2] 
            idx2 = [n2*3, n2*3+1, n2*3+2] 
            
            p1_cur = Coord[:, n1] + U[idx1]
            p2_cur = Coord[:, n2] + U[idx2]
            
            vec_cur = p2_cur - p1_cur
            L_cur = np.linalg.norm(vec_cur)
            
            if L_cur == 0: L_cur = 1e-10 
            cx, cy, cz = vec_cur / L_cur
            
            extension = L_cur - L0_Elem[i]
            equivalent_tension = -dlo[i] * EA_Elem[i] / L0_Elem[i]
            N = (EA_Elem[i] / L0_Elem[i]) * extension * K_ratio[i]
            N += equivalent_tension
            if Elem_Types[i] ==2 or Elem_Types[i] == 4:
                if N<0:
                    N = 0
            equivalent_loads[int(Connect[i][0]*3)] += -equivalent_tension * -cx
            equivalent_loads[int(Connect[i][0]*3+1)] += -equivalent_tension * -cy
            equivalent_loads[int(Connect[i][0]*3+2)] += -equivalent_tension * -cz
            equivalent_loads[int(Connect[i][1]*3)] += -equivalent_tension * cx
            equivalent_loads[int(Connect[i][1]*3+1)] += -equivalent_tension * cy
            equivalent_loads[int(Connect[i][1]*3+2)] += -equivalent_tension * cz

            axial_forces[i] = N
            
            fx, fy, fz = N * cx, N * cy, N * cz
            
            # Assemblage
            P_int[idx1[0]] -= fx; P_int[idx1[1]] -= fy; P_int[idx1[2]] -= fz
            P_int[idx2[0]] += fx; P_int[idx2[1]] += fy; P_int[idx2[2]] += fz
            
            # Masse fictive
            k_val = EA_Elem[i] / L0_Elem[i] * K_ratio[i] # Raideur de l'élément (modifiée par DR)
            for idx in idx1: K_diag_global[idx] += k_val # Simplifié pour la masse
            for idx in idx2: K_diag_global[idx] += k_val

        # Résidu
        #Forces += equivalent_loads
       
        R = Forces - P_int
        R[Fixed_DoF] = 0.0 # Ici Fixed_DoF doit être une liste d'indices !
        
        res_norm = np.linalg.norm(R)
        if res_norm < tol_residual:
            print(f"Convergence à l'etape {step}. Residu : {res_norm:.2e}")
            break
        elif step >= max_steps:
            convergence = False
            break
        M = 2 * (dt**2) * K_diag_global * mass_factor
        M = np.maximum(M, 1e-6) 
        M[Fixed_DoF] = 1e20 # Masse infinie aux appuis
        
        V_inc = dt * (R / M)
        if np.allclose(V, 0, atol=1e-15): V_inc *= 0.5
        V_next = V + V_inc
        
        # Kinetic Damping
        # On calcule l'énergie cinétique uniquement sur les DDL libres
        active_mask = np.ones(No_Ddl, dtype=bool)
        active_mask[Fixed_DoF] = False # Désactive les DDL fixés
        
        KE_next = 0.5 * np.sum(M[active_mask] * V_next[active_mask]**2)
        
        U_correction = np.zeros(No_Ddl)
        if step > 1 and KE_next <= KE_prev:
            U_correction = -1.5 * dt * V_next + 0.5 * (dt**2) * (R / M)
            V_next = 0.5 * dt * R / M
            KE_next = 0.0
            
        U += (dt * V_next) + U_correction
        V = V_next
        KE_prev = KE_next
        KE_list.append(KE_prev)
    print(f"Fin de la simulation apres {step} etapes. Residu final : {res_norm:.2e}")
    
    return U, axial_forces,convergence