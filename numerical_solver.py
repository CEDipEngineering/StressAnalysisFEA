import numpy as np
import time
def jacobi_method(ite, tol, K, F):
    xList = np.zeros(len(K))
    xListNext = xList.copy()
    count = 0
    n = len(K)
    while count < ite:
        # x1 = [F[0] - (K[0][1] * x2 + K[0][2] * x3)] / k[0][0]
        # x2
        for i in range(n):
            xListNext[i] = F[i]
            for j in range(n):
                if j != i:
                    xListNext[i] -= K[i][j] * xList[j]
            xListNext[i] /= K[i][i]
        
        #print((max(abs(xList-xListNext))))
        if(0 not in xList):
            tolCalc = max(abs(np.array([(x-y)/x for x, y in zip(xList, xListNext)])))
            if (tolCalc<tol):
                print(f"Rodei {count} vezes")
                return xList

        xList = xListNext.copy()
        count += 1
    # print(f"Rodei {count} vezes")
    return xList


def gauss_seidel_method(ite, tol, K, F):
    xList = np.zeros(len(K))
    xListNext = xList.copy()
    count = 0
    n = len(K)
    while count < ite:
        #xListNext = xList.copy()
        for i in range(n):
            xListNext[i] = F[i]
            for j in range(n):
                if j != i:
                    xListNext[i] -= K[i][j] * xListNext[j]
            xListNext[i] /= K[i][i]
        
        if(0 not in xList):
            tolCalc = max(abs(np.array([(x-y)/x for x, y in zip(xList, xListNext)])))
            if (tolCalc<tol):
                #print(f"Rodei {count} vezes")
                return xList

        xList = xListNext.copy()
        count += 1
    # print(f"Rodei {count} vezes")
    return xList

if __name__ == "__main__":


    K = [[3,-0.1,-0.2],[0.1,7,-0.3],[0.3,-0.2,10]]        
    F = [7.85, -19.3, 71.4]
    ite = 100
    tol = 1e-15
    print("Jacobi:")
    start = time.perf_counter()
    print(jacobi_method(ite, tol, K, F))
    print(f"Levei {time.perf_counter() - start}")

    print("Gauss-Seidel:")
    start = time.perf_counter()
    print(gauss_seidel_method(ite, tol, K, F))
    print(f"Levei {time.perf_counter() - start}")