import numpy as np

####################################################################################################################################################################

def matrix_inverse(A,b):
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
    #Ainv = np.around(Ainv,5)
    return Ainv,b

################################################################################################################################################################ 

def matrix_inverse_pp(A,b):       
    dim = np.shape(A)[0]
    Ainv = np.identity(dim)
    i = 0
    while i < dim:             #Damit ich gehe ich jede Zeile schrittweise durch
        u = i
        while u < dim-1:                         #Partiell Pivoting
            if A[i][i] < A[u+1][i]:
                A[[i,u+1]] = A[[u+1,i]]
                Ainv[[i,u+1]] = Ainv[[u+1,i]]
                b[[i,u+1]] = b[[u+1,i]]
            u = u + 1

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
    Ainv = np.around(Ainv,5)
    return Ainv,b

################################################################################################################################################################

def lgs_solver(A,b):
    Ainv,b_0 = matrix_inverse_pp(A,b)
    i = 0
    dim = np.shape(A)[0]
    x = np.zeros(dim)

    while i < dim:
        buffer = 0
        j = 0
        while j < dim:
            buffer = buffer + Ainv[i][j]*b_0[j]
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

###################################################################################################################################################################

def matrix_dot(X,Y):

    if len(np.shape(Y)) < 2 and len(np.shape(X)) < 2:               #Vektor mal Vektor
        result = 0
        for i in range(len(X)):
            result += Y[i]*X[i]
    
    elif len(np.shape(Y)) < 2:                                      #Matrix mal Vektor
        result = np.zeros(np.shape(Y))
        for i in range(len(X)):
            for k in range(len(Y)):
                result[i] += X[i][k] * Y[k]
    else:                                                            #Matrix mal Matrix
     
        row = np.shape(X)[0]

        coloumn = np.shape(Y)[1]
        result = np.zeros(row,coloumn)
        for i in range(len(X)):
            for j in range(len(Y[0])):
                for k in range(len(Y)):
                    result[i][j] += X[i][k] * Y[k][j]

    return result

##################################################################################################################################################################

def crout(A):
    dim = np.shape(A)[0]
    L = np.identity(dim)
    R = np.zeros((dim,dim))

    i = 0

    while i < dim:

        j = i
        while j < dim:
            sum = 0
            for k in range(i):
                sum += L[i][k]*R[k][j]
            R[i][j] = A[i][j] - sum
            j += 1

        j = i

        p = i
        while p < dim:
            sum = 0
            for k in range(j):
                sum += L[p][k]*R[k][j]
            L[p][j] = (1/R[j][j])*(A[p][j] - sum)
            p += 1
        i += 1

    return L,R

######################################################################################################################################################################

def lgs_solver_crout(A,b):
    L,R = crout(A)
    N = np.shape(L)[0]
    y = np.zeros(N)
    y[0] = b[0]/L[0][0]
    i = 1
    while i < N:
        sum = 0
        for j in range(i):
            sum += L[i][j]*y[j]
        y[i] = (1/L[i][i])*(b[i]-sum)
        i += 1
    
    x = np.zeros(N)

    x[N-1] = (y[N-1]/(R[N-1][N-1]))

    i = N - 2

    while i >= 0:
        sum = 0
        j = i + 1
        while j < N:
            sum += R[i][j]*x[j]
            j += 1
        x[i] = (1/R[i][i])*(y[i]-sum)
        i -= 1
    return x

########################################################################################################################################################################

def neville_algorithm(x,m,xdata,ydata):
    size = np.shape(xdata)[0]

    loc = 0
    i = 0
    while i < size-1:
        if xdata[i] <= x <= xdata[i+1]:
            loc = i
        i += 1
    
    if x >= xdata[size-1]:
        loc = size - 1
    
    P = np.asfarray([])
    for i in range(size):                               #Erster Iterationsschritt
        P = np.append(P,ydata[i])

    k = 1                                               #k entspricht der Ordnung der Polynome

    while k < m:                                        #Ab hier werden die höheren Polynome bestimmt, angefangen mit linearer Ordnung (k=1)
        i = 0
        P_higher_order = np.zeros(size-k)
        while (i+k) < size:
            P_higher_order[i] = ((x-xdata[i+k])*P[i]+(xdata[i]-x)*P[i+1])/(xdata[i]-xdata[i+k])
            i += 1
        P = P_higher_order
        k += 1
    
    loc = loc - m + 1

    if loc <= 0:
        loc = 0

    return P[loc]

#############################################################################################################################################################################
    