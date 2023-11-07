import numpy as np

####################################################################################################################################################################

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

################################################################################################################################################################ 

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


################################################################################################################################################################

def tridiag_LU_decomp(a,b,c,vecb):
    N = np.shape(b)[0]
    f = c

    e = np.zeros(N)

    d = np.zeros(N)          #ignorier den ersten Wert

    e[0] = b[0]
    d[1] = a[1]/e[0]

    i = 1

    while i < N:                    #BERECHNUNG VON e,d,f. FUNKTIONIERT
        e[i] = b[i] - d[i]*f[i-1]
        if (i+1) == N:
            i = i + 1
            continue
        d[i+1] = a[i+1]/e[i]
        i += 1

    L = np.identity(N)
    R = np.identity(N)

    j = 0

    while j < N:                                      #EINSETZEN IN MATRIZEN L,R FUNKTIONIERT!
        R[j][j] = e[j]
        if j == N-1:
            j += 1
            continue
        L[j+1][j] = d[j+1]
        R[j][j+1] = f[j]
        j += 1
    
    x = np.zeros(N)

    y = np.zeros(N)

    y[0] = vecb[0]

    i = 1
    while i < N:                          #Berechnung von y, FUNKTIONIERT, HABS GECHECKT
        j = 0
        y[i] = vecb[i]
        while j < i:
            y[i] -= L[i][j]*y[j]
            j += 1
        i += 1

    x[N-1] = y[N-1]/(R[N-1][N-1])

    i = N-2                                 #Berechnung von x, FUNKTIONIERT HABS GETESTET
    while i >= 0:
        j = N-1
        x[i] = y[i]
        while j > i:
            x[i] -=  R[i][j]*x[j]
            j -= 1
        x[i] = x[i]/R[i][i]
        i -= 1
    
    return L,R,x

##################################################################################################################################################################
    