import numpy as np
from numerical_solver import *
ite = 1000
tol = 1e-12
K = np.array([[47.7, 0],[0, 17.7]])
F = np.array([[150], [-200]])
print("Gauss-Seidel:")
start = time.perf_counter()
u = gauss_seidel_method(ite, tol, K, F).tolist() #MODELO USADO NO CALCULO FINAL
print(f"Levei {time.perf_counter() - start}\n\n")
print(u)
x = np.linalg.solve(K, F)
print(x)