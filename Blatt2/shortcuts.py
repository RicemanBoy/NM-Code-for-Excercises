import numpy as np

def matrix_inverse(A):
    dim = np.shape(A)[0]
    Ainv = np.identity(dim)
    i = 0
    while i < dim:             #Damit ich gehe ich jede Zeile schrittweise durch

        if A[i][i] != 0:
            c = A[i][i]                 #c ist der Vorfaktor vor dem (i,i)-Eintrag der Matrix
            A[i] = A[i]/c            #Division der i-ten Zeile durch den Vorfaktor an Stelle i,i
            Ainv[i] = Ainv[i]/c       #Selbiges mit der i-ten Zeile in der Einheitsmatrix

            x = 0

            while x < dim:              #Damit gehe ich jede Zeile durch, die ich von der i-ten Zeile subtrahieren muss
                if x == i:
                    x = x + 1
                    continue

                k = A[x][i]           #Koeffizient mit dem ich die Zeile i multiplizieren muss, damit ich Zeile i von i+1 abziehen kann

                j = 0                   #Mit j gehe ich die Spalten alle durch

                while j < dim:          #Damit ich gehe ich jede Spalte j in EINER Zeile durch
                    A[x][j] = A[x][j] - k*A[i][j]
                    Ainv[x][j] = Ainv[x][j] - k*Ainv[i][j]
                    j = j + 1
                x = x + 1

            i = i + 1
            #print(np.matrix(A))
            #print(np.matrix(Ainv))

        else:
            A[[i,i+1]] = A[[i+1,i]]
            Ainv[[i,i+1]] = Ainv[[i+1,i]]
            b[[i,i+1]] = b[[i+1,i]]
    Ainv = np.around(Ainv,5)
    return Ainv
    
def lgs_solver(A,b):
    Ainv = matrix_inverse(A)
    i = 0
    dim = np.shape(A)[0]
    x = np.zeros(dim)

    while i < dim:
        buffer = 0
        j = 0
        while j < dim:
            buffer = buffer + Ainv[i][j]*b[j]
            j = j + 1
        x[i] += buffer
        i += 1
    x = np.around(x,5)
    return x

#def LU_decomp(A):
    