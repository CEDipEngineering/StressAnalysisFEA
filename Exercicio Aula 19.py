class Element():

    def __init__(self, E, A):
        self.nodes = []
        self.conectivity = []
        self.E = E
        self.A = A
        self.tension = 0
        self.deformation = 0
        self.l = 0

    def connect_nodes(self, nodeList):
        # Devia ser sempre 2
        for nod in nodeList:
            self.nodes.append(nod)

        self.l = np.sqrt((self.nodes[0].x-self.nodes[1].x)**2 + (self.nodes[0].y-self.nodes[1].y)**2)

    def calculateS(self):
        self.m = np.array(self.m)
        self.S = (self.E*self.A/self.l)*(np.dot(self.m[:,None], self.m.transpose()[None,:]))/(np.linalg.norm(self.m)**2)
        self.K = np.kron(np.dot(np.array(self.conectivity)[:,None], np.array(self.conectivity).transpose()[None,:]), self.S)

    def calculateDefTens(self):
        localU = np.array([[x.u, y.v] for x,y in zip(self.nodes, self.nodes)])
        localU = localU.reshape((2*len(self.nodes), 1))

        trigMatrix = [-self.cosT(), -self.sinT(), self.cosT(), self.sinT()]

        self.tension = self.E/self.l * (np.dot(trigMatrix, localU))
        self.tension = self.tension[0]
        self.deformation = 1/self.l * (np.dot(trigMatrix, localU))
        self.deformation = self.deformation[0]

    def sinT(self):
        return ((self.nodes[1].y - self.nodes[0].y) / self.l)
        
    def cosT(self):
        return ((self.nodes[1].x - self.nodes[0].x) / self.l)
        
    def __repr__(self) -> str:
        return f"Connects node {self.nodes[0]} to {self.nodes[1]}"

class Node():
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.u = -1
        self.v = -1

    def __add__(self, other):
        return Node(self.x+other.x, self.y+other.y)

    def __sub__(self, other):
        return Node(self.x-other.x, self.y-other.y)

    def __repr__(self):
        return f"x:{self.x} y:{self.y} u:{self.u} v:{self.v}"

    def asList(self):
        return [self.x, self.y]

if __name__ == "__main__":
    
    import numpy as np
    from funcoesTermosol import importa as LeExcel
    import time
    from numerical_solver import *
    
    nn,N,nm,Inc,nc,F,nr,R = LeExcel('entradaAPS.xlsx')
    nodes = [Node(x, y) for x,y in zip(N[0],N[1])]
    elements = []
    conectivity_g = []
    for conection in Inc:
        el = Element(conection[2], conection[3])
        el.connect_nodes([nodes[int(conection[0])-1], nodes[int(conection[1])-1]])
        el.conectivity = [0]*nn
        el.conectivity[int(conection[0])-1] = -1
        el.conectivity[int(conection[1])-1] = 1
        conectivity_g.append(el.conectivity)
        elements.append(el)
    
    conectivity_g = np.array(conectivity_g).transpose()
    M = np.matmul(N, conectivity_g)
    Me = [[x, y] for x,y in zip(M[0],M[1])]
    K_g = []
    for i in range(nm):
        elements[i].m = Me[i]
        elements[i].calculateS()
        K_g.append(elements[i].K)
    K_g = sum(K_g)
    #print(K_g)

    #print(R)

    K_g_restricted = np.delete(K_g, R, axis=0)
    K_g_restricted = np.delete(K_g_restricted, R, axis=1)
    F_restricted = np.delete(F, R, axis=0)
    
    # print(K_g_restricted)

    ite = 1000
    tol = 1e-12

    print("Jacobi:")
    start = time.perf_counter()
    lixo = jacobi_method(ite, tol, K_g_restricted, F_restricted)
    print(f"Levei {time.perf_counter() - start}")

    print("Gauss-Seidel:")
    start = time.perf_counter()
    u = gauss_seidel_method(ite, tol, K_g_restricted, F_restricted).tolist()
    print(f"Levei {time.perf_counter() - start}\n\n")

    x = np.linalg.solve(K_g_restricted, F_restricted)
    # print(u)
    # print(x)

    # WTF
    for i in range(2*nn):
        if i in R.astype(int):
            u.insert(i, 0)
        else:
            u[i] = u[i]
        if i%2 == 0:
            nodes[i//2].u = u[i]
        else:
            nodes[i//2].v = u[i]

    f_f = K_g.dot(u)
    for i in range(len(f_f)):
        if abs(f_f[i]) < 1e-9:
            f_f[i] = 0
    # print(f_f)

    # print(elements[0].conectivity, np.array(u).transpose())
    for seg in elements:
        seg.calculateDefTens()
    
    print(type(seg.tension))






    nodalString = """--------------------------------------
    Vetor de deslocamento nodal [m]:
    """
    forceString = """--------------------------------------
    Forças [N]:
    """
    tensionString = """--------------------------------------
    Tensão em cada elemento [Pa]:
    """
    deformationString = """--------------------------------------
    Deformação longitudinal específica em cada elemento:
    """

    finalString=""

    finalString += nodalString
    for i in u:
        finalString += f"\n{i:.3e}"
    finalString += "\n\n"

    finalString += forceString
    for i in f_f:
        finalString += f"\n{i:.3e}"
    finalString += "\n\n"

    finalString += tensionString
    for i in elements:
        finalString += f"\n{i.tension:.3e}"
    finalString += "\n\n"

    finalString += deformationString
    for i in elements:
        finalString += f"\n{i.deformation:.3e}"
    finalString += "\n\n"

    with open("out.txt", "wb") as out:
        out.write(finalString.encode("utf8"))