import numpy as np

def treillis(n,H,L):
    Coordonates_base = np.zeros((2*n+1, 3))
    Coordonates_base[0] = [0, 0,0]
    for i in range(n):
        Coordonates_base[i+1] = [(i+1)*L, 0,0]
        Coordonates_base[n+i+1] = [i*L+L/2, H, 0]

    a = 1
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.zeros((4*n-1, 6))
    Connections_base[0] = [0, 1, a, 0,0,0]
    Connections_base[n] = [0, n+1, a, 0,0,0]
    Connections_base[2*n] = [1, n+1, a, 0,0,0]
    for i in range(n-1):
        Connections_base[i+1] = [i+1, i+2, a, 0,0,0]
        Connections_base[n+i+1] = [i+1, n+i+2, a, 0,0,0]
        Connections_base[2*n+i+1] = [i+2, n+i+2, a, 0,0,0]
        Connections_base[3*n+i] = [n+i+1, n+2+i, a, 0,0,0]
    
    Pretension = np.zeros(4*n-1)
    # [noeud, axe] 
    Constraints = np.array([
        [0, "x",0],
        [0, "y",0],
        [0,"z",0],
        [n,"y",0],
        [n,"z",0]
        
        
    ])
    P = -10000
    # [N]
    # [noeud, axe, valeur]
    Forces =[]
    for i in range(n):
        Forces.append([i+n+1, "y", P])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def trois_cables():

    Coordonates_base = np.array([
        [0,0,0],
        [0,0,1],
        [-1,0,1],
        [1,0,1]
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,2,0,0,0],
        [1, 2,2,0,0,0],
        [1, 3,2,0,0,0],
    ])
    Pretension = np.array([-0.10,-0.1,-0.1])
    # [noeud, axe] 
    Constraints = np.array([
        [0, "x",0],
        [0, "y",0],
        [0,"z",0], 
        [2,"x",0],
        [2,"y",0],
        [2,"z",0],
        [3,"x",0],
        [3,"y",0],
        [3,"z",0]
        
    ])
    P = 0
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [1, "z", -P],
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def deux_cables():

    Coordonates_base = np.array([
        [0,0,0],
        [1,0,0],
        [2,0,0]
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,2,0,0,0],
        [1, 2,2,0,0,0]
    ])
    Pretension = np.array([-0.01,-0.01])
    # [noeud, axe] 
    Constraints = np.array([
        [0, "x",0],
        [0, "y",0],
        [0,"z",0],
        [2,"x",0],
        [2,"y",0],
        [2,"z",0]
        
    ])
    P = 1000
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [1, "y", -P],
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension
 
def poutres(e0,L0,N,e):
    Coordonates_base = np.array([
        [0,0,0],
        [L0,0,0]
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,1,e0,0,0]
    ])
    Pretension = np.array([0,0])
    # [noeud, axe,valeur] 
    Constraints = np.array([
        [0, "y",0],
        [0,"z",0],
        [0,"x",0],
        [1,"y",0],
        [1,"z",0]
    ])
    
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [1,"x",-N],
        [1,"rz",N *e],
        [0,"rz",-N * e]
    ])
    #[m]
    # [noeud, axe, valeur]
    

    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def contreventement():
    Coordonates_base = np.array([
        [0,0,0],
        [1,0,0],
        [1,0,1],
        [0,0,1]
    ])
    a = 3
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
        [0, 1,a,0,0,0],
        [1, 2,a,0,0,0],
        [2, 3,a,0,0,0],
        [3, 0,a,0,0,0],
        [0, 2,2,0,0,0],
        [1, 3,2,0,0,0],
    ])
    b = -0.1
    Pretension = np.array([0,0,0,0,b,b])
    # [noeud, axe] 
    Constraints = np.array([
        [0, "x",0],
        [0, "y",0],
        [0,"z",0],
        [1,"x",0],
        [1,"y",0],
        [1,"z",0]
        
    ])
    P = 100000
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [3, "y", -P]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def simplex(P,e0):
    Coordonates_base = np.array([[0.0, -2.0438, 0.0], 
                                 [0.0, 0.0, 0.0], 
                                 [1.77, -1.0219, 0.0], 
                                 [0.590, -2.2019, 1.950], 
                                 [-0.4319, -0.4319, 1.950], 
                                 [1.6119, -0.4319, 1.950]])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), rolule_start, rotule_end]
    a =1
   
    Connections_base = np.array([[3, 0,2,0,0,0],
                                 [5, 2,2,0,0,0],
                                  [4, 1,2,0,0,0],
                                  [0, 1,2,0,0,0],
                                  [1, 2,2,0,0,0],
                                  [2, 0,2,0,0,0],
                                  [3, 4,2,0,0,0],
                                  [4, 5,2,0,0,0],
                                  [5, 3,2,0,0,0],
                                  [5, 0,a,100,0,0], 
                                  [3, 1,a,e0,0,0], 
                               [4, 2,a,e0,0,0]])
    b = -0.00122
    Pretension = np.array([b,b,b,0,0,0,0,0,0,0,0,0])
    # [noeud, axe] 
    Constraints = np.array([
        [0, "y",0],
        [0, "z",0],
        [1,"x",0],
        [1,"y",0],
        [1,"z",0],
        [2,"z",0]
    ])
    
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [3, "z", -P],
        [4, "z", -P],
        [5, "z", -P]

    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def tensegrité2_2D():
   
    Coordonates_base = np.array([
        [0,0,0],
        [1,-0.2,0],
        [0.8,0.8,0],
        [-0.2,1,0],
        [2.2,1,0],
        [2,0,0]
    ])
    
    a = 3
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    b = -0.1
    e0 = 0
    Connections_base = np.array([
    
        [0, 2,a,e0,1,1],
        [3, 1,a,e0,1,1],
        [1, 4,a,e0,1,1],
        [2, 5,a,e0,1,1],
        [0, 1,2,0,0,0],
        [0, 3,2,0,0,0],
        [1,2,2,0,0,0],
        [2,3,2,0,0,0],
        [2,4,2,0,0,0],
        [1,5,2,0,0,0],
        [4,5,2,0,0,0]])
    Pretension = np.array([0,0,0,0,0,0,-0.06,0,0,0,0])
    # [noeud, axe] 
    Constraints = np.array([
        [0, "x",0],
        [0, "y",0],
        [0,"z",0],
        
        [5,"y",0],
        [5,"z",0]
        
    ])
    P = 70000
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [2, "y", -P]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def poutres_snap_through(P):
    Coordonates_base = np.array([
        [-1,0,0],
        [0,0,0.1],
        [1,0,0]
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,1,0,0,0],
        [1,2,1,0,0,0]

    ])
    Pretension = np.array([0,0])
    # [noeud, axe] 
    Constraints = np.array([
        [0,"x",0],
        [0, "y",0],
        [0,"z",0],
        [0,"rx",0],
        [2,"x",0],
        [2,"z",0],
        [2,"y",0]
        
    ])
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [1,"z",-30000],
        [2,"x",-P],
        [0,"x",P]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def passerelle_deux_modules(P):
    Coordonates_base = np.array([
        [0,-2.163,0.766],
        [0,-0.418,2.256],
        [0,0,0],
        [2.944,-1.721,0],
        [2.944,-1.303,2.256],
        [2.944,0.442,0.766],
        [5.888,-2.163,0.766],
        [5.888,-0.418,2.256],
        [5.888,0,0]
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,4,0,0,0],
        [1,2,4,0,0,0],
        [0,2,4,0,0,0],
        [0,3,4,0,0,0],
        [1,4,4,0,0,0],
        [2,5,4,0,0,0],
        [0,5,3,0,1,1],
        [1,3,3,0,1,1],
        [2,4,3,0,1,1],
        [3,4,4,0,0,0],
        [4,5,4,0,0,0],
        [5,3,4,0,0,0],
        [3,6,4,0,0,0],
        [4,7,4,0,0,0],
        [5,8,4,0,0,0],
        [3,7,3,0,1,1],
        [4,8,3,0,1,1],
        [5,6,3,0,1,1],
        [6,7,4,0,0,0],
        [7,8,4,0,0,0],
        [8,6,4,0,0,0]
        
        
            ])
    a = -0.051
    Pretension = np.array([0,0,0,
                           a,0,a,
                           0,0,0,
                           0,0,0,
                           a,0,a,
                           0,0,0,
                           0,0,0,])
    # [noeud, axe] 
    Constraints = np.array([
        [0,"z",0],
        [2,"x",0],
        [2, "y",0],
        [2,"z",0],
        [6,"z",0],
        [8,"y",0],
        [8,"z",0]
        
    ])
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [3,"z",-P],
        [5,"z",-P]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def poutre_cantilever(M):
    Coordonates_base = np.array([
        [0,0,0],
        [10,0,0]
        
    ])
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
    
        [0, 1,1,0,0,0]

    ])
    Pretension = np.array([0])
    # [noeud, axe] 
    Constraints = np.array([
        [0,"x",0],
        [0, "y",0],
        [0,"z",0],
        [0,"rx",0],
        [0,"ry",0],
        [0,"rz",0]
        
    ])
    
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [1,"ry",M]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

def passerelle_cinq_modules(P,e0):
    Coordonates_base = np.array([
        [0,-2.163,0.766], #0
        [0,-0.418,2.256], #1
        [0,0,0], #2
        [2.944,-1.721,0], #3
        [2.944,-1.303,2.256], #4
        [2.944,0.442,0.766], #5
        [5.888,-2.163,0.766], #6
        [5.888,-0.418,2.256], #7
        [5.888,0,0], #8
        [8.832,-1.721,0], #9
        [8.832,-1.303,2.256], #10
        [8.832,0.442,0.766], #11
        [11.776,-2.163,0.766], #12
        [11.776,-0.418,2.256], #13
        [11.776,0,0], #14
        [14.720,-1.721,0], #15
        [14.720,-1.303,2.256], #16
        [14.720,0.442,0.766] #17

    ])
    r = 1
    # [noeud1, noeud2, type_elem, imperfection(L0/e0,0 si pas d'imperfection), prétension [N]]
    Connections_base = np.array([
        [2,3,4,0,0,0],  #câble additionnel
        [0, 1,4,0,0,0],
        [1,2,4,0,0,0],
        [0,2,4,0,0,0],
        [0,3,4,0,0,0],
        [1,4,4,0,0,0],
        [2,5,4,0,0,0],
        [0,5,3,e0,r,r],
        [1,3,3,e0,r,r],
        [2,4,3,e0,r,r],
        [3,8,4,0,0,0], #câble additionnel
        [3,4,4,0,0,0],
        [4,5,4,0,0,0],
        [5,3,4,0,0,0],
        [3,6,4,0,0,0],
        [4,7,4,0,0,0],
        [5,8,4,0,0,0],
        [3,7,3,e0,r,r],
        [4,8,3,e0,r,r],
        [5,6,3,e0,r,r],
        [8,9,4,0,0,0], #câble additionnel
        [6,7,4,0,0,0],
        [7,8,4,0,0,0],
        [8,6,4,0,0,0],
        [6,9,4,0,0,0],
        [7,10,4,0,0,0],
        [8,11,4,0,0,0],
        [6,11,3,e0,r,r],
        [7,9,3,e0,r,r],
        [8,10,3,e0,r,r],
        [9,14,4,0,0,0], #câble additionnel
        [9,10,4,0,0,0],
        [10,11,4,0,0,0],
        [11,9,4,0,0,0],
        [9,12,4,0,0,0],
        [10,13,4,0,0,0],
        [11,14,4,0,0,0],
        [9,13,3,e0,r,r],
        [10,14,3,e0,r,r],
        [11,12,3,e0,r,r],
        [14,15,4,0,0,0], #câble additionnel
        [12,13,4,0,0,0],
        [13,14,4,0,0,0],
        [14,12,4,0,0,0],
        [15,12,4,0,0,0],
        [16,13,4,0,0,0],
        [17,14,4,0,0,0],
        [12,17,3,e0,r,r],
        [13,15,3,e0,r,r],
        [14,16,3,e0,r,r],
        [15,16,4,0,0,0],
        [16,17,4,0,0,0],
        [17,15,4,0,0,0]
        
        
            ])
    a = -0.028
    b = -0.0533
    c = -0.0693
    Pretension = np.array([0,0,0,0,
                           a,0,a,
                           0,0,0,
                           0,0,0,0,
                           b,0,b,
                           0,0,0,
                           0,0,0,0,
                           c,0,c,
                           0,0,0,
                           0,0,0,0,
                           b,0,b,
                           0,0,0,
                           0,0,0,0,
                           a,0,a,
                           0,0,0,
                           0,0,0
                           ])
    # [noeud, axe] 
    Constraints = np.array([
        
        [0,"z",0],
        [2,"x",0],
        [2, "y",0],
        [2,"z",0],
        [17,"z",0],
        [15,"y",0],
        [15,"z",0]
        
    ])
    # [N]
    # [noeud, axe, valeur]
    Forces = np.array([
        [7,"z",-P],
        [10,"z",-P]
    ])
    return Coordonates_base, Connections_base, Constraints, Forces, Pretension

