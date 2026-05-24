import numpy as np
import matplotlib.pyplot as plt
import math

def create_mesh(init_mesh, Coordonates_base, Connections_base, Constraints, Forces, pretensions, Elem_properties, PP):
    points_dict = {}
    x_coord = []
    y_coord = []
    z_coord = []
    connections = []
    Elem_Types = []
    rotule_start = [] 
    rotule_end = []  
    Pretension = []
    
    # Conversion des contraintes de base en tenant compte du mesh
    Constraints_2 = np.zeros(len(Coordonates_base) * 6)
    Forces_2 = np.zeros(len(Coordonates_base) * 6)
    Deplacements = np.zeros(len(Coordonates_base) * 6)
    
    for i in range(len(Constraints)):
        index = int(Constraints[i][0]) * 6
        if Constraints[i][1] == "x":     
            index += 0
        elif Constraints[i][1] == "y":
            index += 1
        elif Constraints[i][1] == "z":
            index += 2
        elif Constraints[i][1] == "rx":
            index += 3
        elif Constraints[i][1] == "ry":
            index += 4
        elif Constraints[i][1] == "rz":
            index += 5
        Constraints_2[index] = 1
        Deplacements[index] = Constraints[i][2]

    # Conversion des forces de base en tenant compte du mesh
    for i in range(len(Forces)):
        index = int(Forces[i][0]) * 6
        if Forces[i][1] == "x":     
            index += 0
        elif Forces[i][1] == "y":
            index += 1
        elif Forces[i][1] == "z":
            index += 2
        elif Forces[i][1] == "rx":
            index += 3 
        elif Forces[i][1] == "ry":
            index += 4
        elif Forces[i][1] == "rz":
            index += 5
        Forces_2[index] = Forces[i][2]
  
    Constraints_3 = []
    Forces_3 = []
    Initial_point_index = []
    U = []
    
    for i in range(len(Connections_base)): 
        conn = Connections_base[i]
        type_elem = Connections_base[i][2] # Renommé 'type_elem' pour éviter de surcharger le mot-clé 'type'
        prétension = pretensions[i]
        start = Coordonates_base[int(conn[0])]
        end = Coordonates_base[int(conn[1])]

        # -- MODIFICATION ICI --
        # Sortir les calculs globaux de l'élément (indépendants du maillage)
        alpha = np.arctan2(end[1] - start[1], end[0] - start[0])
        theta = np.arctan2(end[2] - start[2], np.sqrt((end[0] - start[0])**2 + (end[1] - start[1])**2))
        L0 = np.sqrt((end[0] - start[0])**2 + (end[1] - start[1])**2 + (end[2] - start[2])**2)

        # Création des segments (start_t, end_t, n_mesh, type_du_segment)
        segments = []
        
        # Hypothèse: les coordonnées sont en mètres, donc 14 cm = 0.14. 
        # (Si tu es en mm, remplace 0.14 par 140 et 0.28 par 280)
        if type_elem == 3 and L0 > 0.28: 
            t_14 = 0.14 / L0
            segments.append((0.0, t_14, 3, 5))                # Premier bout: 14cm, 3 éléments, type 5
            segments.append((t_14, 1.0 - t_14, init_mesh, 3)) # Milieu: le reste, init_mesh éléments, type 3
            segments.append((1.0 - t_14, 1.0, 3, 5))          # Dernier bout: 14cm, 3 éléments, type 5
        else:
            n_mesh = init_mesh if (type_elem == 1 or type_elem == 3 or type_elem ==6) else 1
            segments.append((0.0, 1.0, n_mesh, type_elem))
            
        total_sub_elems = sum(seg[2] for seg in segments)
        current_sub_elem = 0

        # Boucle sur les différents segments d'un même élément
        for seg_start, seg_end, seg_mesh, current_type in segments:
            for j in range(seg_mesh): # Utilisation de 'j' pour ne pas écraser le 'i' de la boucle principale
                
                # Paramétrisation t0 et t1 à l'échelle globale de l'élément (0 à 1)
                t0 = seg_start + (j / seg_mesh) * (seg_end - seg_start)
                t1 = seg_start + ((j + 1) / seg_mesh) * (seg_end - seg_start)

                if conn[3] == 0:
                    e0 = 0
                else:
                    e0 = L0 / conn[3]
                
                # Coordonnées locales
                x0 = round(start[0] + t0 * L0 * np.cos(alpha) * np.cos(theta) + e0 * np.sin(np.pi * t0)*np.sin(alpha), 6)
                x1 = round(start[0] + t1 * L0 * np.cos(alpha) * np.cos(theta) + e0 * np.sin(np.pi * t1)*np.sin(alpha), 6)
                y0 = round(start[1] + t0 * L0 * np.sin(alpha) * np.cos(theta)- e0 * np.sin(np.pi * t0)*np.cos(alpha), 6)
                y1 = round(start[1] + t1 * L0 * np.sin(alpha) * np.cos(theta)- e0 * np.sin(np.pi * t1)*np.cos(alpha), 6)
                z0 = round(start[2] + t0 * L0 * np.sin(theta), 6)
                z1 = round(start[2] + t1 * L0 * np.sin(theta), 6)

                tmp = add_point(x0, y0, z0, points_dict, x_coord, y_coord, z_coord)
                idx0 = tmp[0]

                if abs(x0 - start[0]) < 1e-8 and abs(y0 - start[1]) < 1e-8 and abs(z0 - start[2]) < 1e-8 and tmp[1] == True:
                    Constraints_3.extend(Constraints_2[int(conn[0]) * 6 + k] for k in range(6))
                    U.extend(Deplacements[int(conn[0]) * 6 + k] for k in range(6))
                    Forces_3.extend(Forces_2[int(conn[0]) * 6 + k] for k in range(6))
                    Initial_point_index.append(idx0)
                elif tmp[1] == True:  
                    Constraints_3.extend(0 for _ in range(6))
                    Forces_3.extend(0 for _ in range(6))
                    U.extend(0 for _ in range(6))

                tmp = add_point(x1, y1, z1, points_dict, x_coord, y_coord, z_coord)
                idx1 = tmp[0]

                if abs(x1 - end[0]) < 1e-6 and abs(y1 - end[1]) < 1e-6 and abs(z1 - end[2]) < 1e-6 and tmp[1] == True:
                    Constraints_3.extend(Constraints_2[int(conn[1]) * 6 + k] for k in range(6))
                    Forces_3.extend(Forces_2[int(conn[1]) * 6 + k] for k in range(6))
                    U.extend(Deplacements[int(conn[1]) * 6 + k] for k in range(6))
                    Initial_point_index.append(idx1)
                elif tmp[1] == True:
                    Constraints_3.extend(0 for _ in range(6))
                    Forces_3.extend(0 for _ in range(6))
                    U.extend(0 for _ in range(6))

                connections.append([idx0, idx1])
                Elem_Types.append(current_type) # Ajout du type spécifique (3 ou 5)
                Pretension.append(prétension)
                
                # Gestion des rotules sur la base du nombre total de sous-éléments
                if current_sub_elem == 0:
                    rotule_start.append(conn[4])
                    rotule_end.append(0)
                elif current_sub_elem == total_sub_elems - 1:
                    rotule_end.append(conn[5])
                    rotule_start.append(0)
                else:
                    rotule_start.append(0)
                    rotule_end.append(0)
                
                current_sub_elem += 1

    # ajout du poids propre au point
    for i in range(len(connections)):
        idx_start = connections[i][0]
        idx_end = connections[i][1]
        length = np.sqrt((x_coord[idx_end]-x_coord[idx_start])**2 + (y_coord[idx_end]-y_coord[idx_start])**2 + (z_coord[idx_end]-z_coord[idx_start])**2)
        surface = Elem_properties[int(Elem_Types[i])-1][0]        
        density = Elem_properties[int(Elem_Types[i])-1][5] 
        
        weight = length * surface * density * 9.81 * PP
        Forces_3[idx_start * 6 + 2] += -weight / 2
        Forces_3[idx_end * 6 + 2] += -weight / 2 
        
    Coord = np.array([x_coord, y_coord, z_coord]) 
    Connect = np.array(connections)
    connections.clear()
     
    Fixed_Dof =[]
    for i in range(len(Constraints_3)):
        if Constraints_3[i] == 1:
            Fixed_Dof.append(i)
    Fixed_Dof = np.array(Fixed_Dof)
    Forces_3 = np.array(Forces_3)
    U = np.array(U)
    
    return Coord, Connect, Fixed_Dof, Forces_3, Elem_Types, Pretension, Initial_point_index, rotule_start, rotule_end, U

def add_point(x, y, z, points_dict, x_coord, y_coord, z_coord):
    key = (x,y,z)
    New_point = False
    if key not in points_dict:
        points_dict[key] = len(x_coord)
        x_coord.append(x)
        y_coord.append(y)
        z_coord.append(z)
        New_point = True
    return points_dict[key], New_point