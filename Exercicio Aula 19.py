from typing import final
from funcoesTermosol import plota


class Element():

    def __init__(self, E, A):
        self.nodes = []
        self.conectivity = []
        self.E = E  #MODULO DE ELASTICIDADE
        self.A = A  #AREA DE SECCAO
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
        self.K = np.kron(np.dot(np.array(self.conectivity)[:,None], np.array(self.conectivity).transpose()[None,:]), self.S) #MATRIZ DE RIGIDEZ DO ELEMENTO

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
    
    nn,N,nm,Inc,nc,F,nr,R = LeExcel('entradasimples.xlsx')
    F*=1

    
    #######################################################################
    # nn = número de nós                                                  #
    # N = Lista dos nós [(x1,y1),(x2,y2),etc] shape = (2,nn)              #
    # nm = número de membros                                              #
    # Inc = Matriz de incidência (qual nó conecta em qual) shape = (2,nm) #
    # nc = número de cargas                                               #
    # F = Matriz de Forças ([nó1, x/y, carga], etc) shape = (3,nc)        #
    # nr = número de restrições                                           #
    # R = Matriz de restrições ([nó1, x/y], etc) shape = (2, nr)          #
    #######################################################################


    density = 848
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


    K_g_restricted = np.delete(K_g, R.astype(int), axis=0)
    K_g_restricted = np.delete(K_g_restricted, R.astype(int), axis=1)
    F_restricted = np.delete(F, R.astype(int), axis=0)
    
    #M = Matriz dos Membros (matriz de nós x matriz de conectividade)
    #Me = Matriz dos Membros de um elemento específico
    #C =  Matriz de Conectividade
    #K = Matriz de Rigidez do elemento
    #Kg = Soma de todos os Ks
    #u = Vetor de deslocamento nodal

    ite = 1000
    tol = 1e-12

    # print("Jacobi:")
    # start = time.perf_counter()
    # lixo = jacobi_method(ite, tol, K_g_restricted, F_restricted)
    # print(f"Levei {time.perf_counter() - start}")

    # print("Gauss-Seidel:")
    # start = time.perf_counter()
    u = gauss_seidel_method(ite, tol, K_g_restricted, F_restricted).tolist() #MODELO USADO NO CALCULO FINAL
    # print(f"Levei {time.perf_counter() - start}\n\n")

    # u = np.linalg.solve(K_g_restricted, F_restricted)


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

    for seg in elements:
        seg.calculateDefTens()
    
    # print(type(seg.tension))



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

    #print(len(elements))

    displacementList = []
    tensionList = []
    deformationList = []
    lengthList = [i.l for i in elements]
    mass = sum([i.l*i.A*density for i in elements])

    finalString += nodalString
    i = 0
    while i < len(u):
        x = f"{u[i]:.3e}"
        y = f"{u[i+1]:.3e}"
        finalString += "\n(" + x + " , " + y + ")"
        displacementList.append(u[i])
        displacementList.append(u[i+1])
        i+=2
    finalString += "\n\n"

    finalString += forceString
    i = 0
    while i < len(f_f):
        x =  f"{f_f[i]:.3e}"
        y =  f"{f_f[i+1]:.3e}"
        finalString += "\n(" + x + " , " + y + ")"
        i += 2
    finalString += "\n\n"

    finalString += tensionString
    i = 0
    while i < len(elements):
        # x =  f"{elements[i].tension:.3e}"
        # y =  f"{elements[i+1].tension:.3e}"
        finalString += f"\n{elements[i].tension:.3e}"
        tensionList.append(elements[i].tension)
        i += 1
    finalString += "\n\n"

    finalString += deformationString
    i = 0
    while i < len(elements):
        # x =  f"{elements[i].deformation:.3e}"
        # y =  f"{elements[i+1].deformation:.3e}"
        finalString += f"\n{elements[i].deformation:.3e}"
        deformationList.append(elements[i].deformation)
        i += 1
    finalString += "\n\n"

    finalString += f"""=======================\nSummary:\n===========================\n\n"""
    finalString += f"Max:\n\tTension: {max(tensionList):.3e}\n\tDeformation: {max(deformationList):.3e}\n\tDisplacement: {max(displacementList):.3e}\n"
    finalString += f"Min:\n\tTension: {min(tensionList):.3e}\n\tDeformation: {min(deformationList):.3e}\n\tDisplacement: {min(displacementList):.3e}\n"
    finalString += f"Força total: {abs(int(sum(F)))}[N]\n"
    finalString += f"Maior comprimento: {abs(max(lengthList))}[m]\n"
    finalString += f"Massa total: {mass*1000}[g]\n"
    finalString += "\n"
    failed = False
    if (max(tensionList)>18e6 or min(tensionList)<-18e6):
        failed = True
        finalString+="Limite de ruptura excedido!\n"
    elif (max(displacementList)>0.02 or min(displacementList)<-0.02):
        failed = True
        finalString+="Limite de deslocamento atingido!\n"
    elif (max(deformationList)>5e-2):
        failed = True
        finalString+="Limite de deformação atingido!\n"
    elif (mass>0.25):
        failed = True
        finalString+="Limite de massa atingido!\n"
    elif abs(max(lengthList))>0.11:
        failed = True
        finalString+="Limite de comprimento atingido!\n"
        
    if failed:
        finalString += "\nA estrutura falhou!\n"
    else:
        finalString += "\nA estrutura está de pé!\n"
    
    with open("out.txt", "wb") as out:
        out.write(finalString.encode("utf8"))


    plota(N, Inc)

    ###############################################################################################
    #Interpretação Out.txt
    #Vetor de deslocamento nodal [m]:Quanto cada nó se deslocou (x,y) em relação ao início
    #Forças [N]: Força de reação em cada nó (x,y)
    #Tensão em cada elemento [Pa]: (x,y)
    #Deformação longitudinal específica em cada elemento: "porcentagem" de mudança de tamanho (x,y)
