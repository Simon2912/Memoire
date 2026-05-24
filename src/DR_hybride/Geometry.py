import numpy as np
import math

def create_mesh_truss(Coordonates_base,Connections_base, Constraints, Forces):
    
    # 1. Création du masque (tableau de 0 et 1)
    Fixed_Dof_Mask = np.zeros(len(Coordonates_base) * 3, dtype=int)
    Forces_Vector = np.zeros(len(Coordonates_base) * 3, dtype=float)    
    
    for i in range(len(Constraints)):
        index = int(Constraints[i][0]) * 3
        if Constraints[i][1] == "x": index += 0
        elif Constraints[i][1] == "y": index += 1
        elif Constraints[i][1] == "z": index += 2
        
        Fixed_Dof_Mask[index] = 1 
    
    # CORRECTION CRITIQUE : Conversion du masque en liste d'indices
    # np.where renvoie un tuple, on prend le premier élément [0]
    Fixed_Dof_Indices = np.where(Fixed_Dof_Mask == 1)[0]

    for i in range(len(Forces)):
        index = int(Forces[i][0]) * 3
        if Forces[i][1] == "x": index += 0
        elif Forces[i][1] == "y": index += 1
        elif Forces[i][1] == "z": index += 2
        Forces_Vector[index] = float(Forces[i][2])
        
    # Elem_types fictif pour l'exemple (tout type 1)
    Elem_types = np.zeros(len(Connections_base), dtype=int)
    initial_point_index = []
    for i in range(len(Connections_base)): # Car 3 éléments dans trois_cables()
        Elem_types[i] = Connections_base[i][2]
    for i in range(len(Coordonates_base)):
        initial_point_index.append(int(i))
    
    return Fixed_Dof_Indices, Forces_Vector, Elem_types, initial_point_index




def create_mesh_beam(init_mesh, Coordonates_base, Connections_base, Constraints, Forces, pretensions,Elem_properties,PP,beta):
    beta = math.radians(beta)  # Convertir l'angle en radians pour les calculs trigonométriques
    points_dict = {}
    x_coord = []
    y_coord = []
    z_coord = []
    connections = []
    Elem_Types = []# Liste pour stocker les types d'éléments, si nécessaire
  
    Pretension = []  # Liste pour stocker les prétensions, si nécessaire

    # Conversion des contraintes de base en tenant compte du mesh
    Constraints_2 = np.zeros(len(Coordonates_base) * 6, dtype=int)  # Initialiser un tableau de zéros
    Forces_2 = np.zeros(len(Coordonates_base) * 6, dtype=int)  # Initialiser un tableau de zéros pour les forces    
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
        Constraints_2[index] = 1  # Marquer la contrainte comme active
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
        Forces_2[index] = Forces[i][2]  # Appliquer la force
  

    
    Constraints_3 = []
    Forces_3 = []
    Initial_point_index = []
    
    for i in range(len(Connections_base)): 

        conn = Connections_base[i]
        type = Connections_base[i][2]
        prétension = pretensions[i]
        start = Coordonates_base[int(conn[0])]
        end = Coordonates_base[int(conn[1])]
        if type == 1 or type == 3:
            n_mesh = init_mesh
        elif type == 2 or type == 4:
            n_mesh = 1
        for i in range(n_mesh):
            t0 = i / n_mesh
            t1 = (i + 1) / n_mesh

            # Angle dans le plan XY
            alpha = np.arctan2(end[1] - start[1], end[0] - start[0])
            # Angle d'inclinaison dans le plan vertical
            theta = np.arctan2(end[2] - start[2], np.sqrt((end[0] - start[0])**2 + (end[1] - start[1])**2))

            # Longueur totale entre les deux points
            L0 = np.sqrt((end[0] - start[0])**2 + (end[1] - start[1])**2 + (end[2] - start[2])**2)
            if conn[3] == 0:
                e0 = 0
            else:
                e0 = L0 / conn[3]
            # Coordonnées locales selon la progression le long de l’élément
            x0 = round(start[0] + t0 * L0 * np.cos(alpha) * np.cos(theta) + e0 * np.sin(np.pi * t0)*np.sin(alpha), 6)
            x1 = round(start[0] + t1 * L0 * np.cos(alpha) * np.cos(theta) + e0 * np.sin(np.pi * t1)*np.sin(alpha), 6)

            y0 = round(start[1] + t0 * L0 * np.sin(alpha) * np.cos(theta), 6)
            y1 = round(start[1] + t1 * L0 * np.sin(alpha) * np.cos(theta), 6)

            z0 = round(start[2] + t0 * L0 * np.sin(theta)- e0 * np.sin(np.pi * t0)*np.cos(alpha), 6)
            z1 = round(start[2] + t1 * L0 * np.sin(theta)- e0 * np.sin(np.pi * t1)*np.cos(alpha), 6)

            # Ajout des points dans le maillage (fonction utilisateur)
            tmp = add_point(x0, y0, z0, points_dict, x_coord, y_coord, z_coord)
            


            idx0 = tmp[0]

            if abs(x0 - start[0]) < 1e-8 and abs(y0 - start[1]) < 1e-8 and abs(z0 - start[2]) < 1e-8 and tmp[1] == True:
                Constraints_3.extend(Constraints_2[int(conn[0]) * 6 + i] for i in range(6))
                Forces_3.extend(Forces_2[int(conn[0]) * 6 + i] for i in range(6))
                Initial_point_index.append(idx0)

            elif tmp[1] == True:  # Si le point est nouveau, on ajoute des zéros
                
                Constraints_3.extend(0 for _ in range(6))
                Forces_3.extend(0 for _ in range(6))

            tmp = add_point(x1, y1, z1, points_dict, x_coord, y_coord, z_coord)
            idx1 = tmp[0]

            if abs(x1 - end[0]) < 1e-6 and abs(y1 - end[1]) < 1e-6 and abs(z1 - end[2]) < 1e-6 and tmp[1] ==True:
                Constraints_3.extend(Constraints_2[int(conn[1]) * 6 + i] for i in range(6))
                Forces_3.extend(Forces_2[int(conn[1]) * 6 + i] for i in range(6))
                Initial_point_index.append(idx1)

            elif tmp[1] == True:
                Constraints_3.extend(0 for _ in range(6))
                Forces_3.extend(0 for _ in range(6))

            connections.append([idx0, idx1])
            Elem_Types.append(type)
            Pretension.append(prétension)
       
    #ajout du poids propre au point
    for i in range(len(connections)):
        idx_start = connections[i][0]
        idx_end = connections[i][1]
        length = np.sqrt((x_coord[idx_end]-x_coord[idx_start])**2 + (y_coord[idx_end]-y_coord[idx_start])**2 + (z_coord[idx_end]-z_coord[idx_start])**2)
        surface = Elem_properties[int(Elem_Types[i])-1][0]  # Récupérer la surface de la section transversale        
        density = Elem_properties[int(Elem_Types[i])-1][5]  # Récupérer la densité
        
        weight = length * surface * density * 9.81* PP  # Calculer le poids propre
        
        Forces_3[idx_start * 6 + 2] += -weight / 2 *math.cos(beta)  # Répartir la moitié du poids au point de départ
        Forces_3[idx_end * 6 + 2] += -weight / 2  *math.cos(beta)  # Répartir l'autre moitié au point d'arrivée
        Forces_3[idx_start * 6 ] += -weight / 2 *math.sin(beta)  # Répartir la moitié du poids au point de départ
        Forces_3[idx_end * 6 ] += -weight / 2  *math.sin(beta)  # Répartir l'autre moitié au point d'arrivée
     # Conversion des listes en tableaux numpy


    Coord = np.array([x_coord, y_coord, z_coord])  # 3 x N
    # Conversion en tableau numpy
    

    Connect = np.array(connections)
    connections.clear()
    

    #Elem_Types = np.array(Elem_Types)  # Convertir en tableau numpy si nécessaire
 
    Fixed_Dof =[]
    for i in range(len(Constraints_3)):
        if Constraints_3[i] == 1:
            Fixed_Dof.append(i)
    Fixed_Dof = np.array(Fixed_Dof)  # Convertir en tableau numpy
    Forces_3 = np.array(Forces_3)
    

      # Convertir en tableau numpy
    return Coord, Connect, Fixed_Dof, Forces_3, Elem_Types, Pretension, Initial_point_index


# Fonction pour ajouter un point s'il n'existe pas encore
def add_point(x, y, z,points_dict, x_coord, y_coord, z_coord):
    key = (x,y,z)  # éviter les problèmes d'arrondi
    New_point = False
    if key not in points_dict:
        points_dict[key] = len(x_coord)
        x_coord.append(x)  # Ajouter x à la première ligne
        y_coord.append(y)
        z_coord.append(z)
        New_point = True
    return points_dict[key],New_point
