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
        result = np.zeros((row,coloumn))
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
    m += 1
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
    
    loc = loc - m + 2

    if loc <= 0:
        loc = 0
        
    if loc > np.size(P)-1:
        loc = np.size(P)-1

    return P[loc]

#############################################################################################################################################################################
    

def cubic_spline_interpolation(x,xdata,ydata):
    N = np.shape(xdata)[0]                             #Anzahl der Stützpunkte xdata_i
    
    loc = 0

    i = 0

    while i < N-1:
        if xdata[i] <= x < xdata[i+1]:
            loc = i
        i += 1

    if x >= xdata[N-1]:
        loc = N - 2                                      #loc ist bei unterem Stützpunkt, also bei j

    
    matrix = np.zeros((N-2,N-2))
    b = np.zeros(N-2)

    #Fall: j = 1, erste Zeile der Matrix mit nur 2 Elementen, statt 3

    matrix[0][0] = (xdata[2]-xdata[0])/3                            #y1''
    matrix[0][1] = (xdata[2]-xdata[1])/6                            #y2''
    b[0] = ((ydata[2]-ydata[1])/(xdata[2]-xdata[1]))-((ydata[1]-ydata[0])/(xdata[1]-xdata[0]))

    j = 2

    i = 1
    while i < N-3:
        matrix[i][i-1] = (xdata[j]-xdata[j-1])/6
        matrix[i][i] = (xdata[j+1]-xdata[j-1])/3
        matrix[i][i+1] = (xdata[j+1]-xdata[j])/6
        b[i] = ((ydata[j+1]-ydata[j])/(xdata[j+1]-xdata[j]))-((ydata[j]-ydata[j-1])/(xdata[j]-xdata[j-1]))
        j += 1
        i += 1

    
    #Fall: j = N - 3, letzte Zeile der Matrix mit nur 2 Elementen, statt 3

    matrix[N-3][N-4] = (xdata[N-2]-xdata[N-3])/6
    matrix[N-3][N-3] = (xdata[N-1]-xdata[N-3])/3
    b[N-3] = ((ydata[N-1]-ydata[N-2])/(xdata[N-1]-xdata[N-2]))-((ydata[N-2]-ydata[N-3])/(xdata[N-2]-xdata[N-3]))
 
    
    y2 = np.asfarray([0])
    y2 = np.append(y2, lgs_solver_crout(matrix,b))
    y2 = np.append(y2,0)

    if loc > np.size(y2)-2:
        loc = np.size(y2) - 2

    A = (xdata[loc+1]-x)/(xdata[loc+1]-xdata[loc])
    B = 1 - A

    ylinear = A*ydata[loc] + B*ydata[loc+1]

    C = (1/6)*((A**3)-A)*((xdata[loc+1]-xdata[loc])**2)
    D = (1/6)*((B**3)-B)*((xdata[loc+1]-xdata[loc])**2)

    p = C*y2[loc] + D*y2[loc+1]

    y = ylinear + p
    
    return y

#######################################################################################################################################################################################

def integral_ext_trap(f,a,b,N):                     #funktioniert
    if b >= a:
        h = (b-a)/N
        x = a
        area = 0.5*(f(a)+f(b))
        i = 1
        while i < N:
            area += f(x+i*h)
            i += 1
        area = h*area
    else:
        h = (a-b)/N
        x = b
        area = 0.5*(f(a)+f(b))
        i = 1
        while i <= N:
            area += f(x+i*h)
            i += 1
        area = -h*area
    return area

#######################################################################################################################################################################################

def integral_ext_simp(f,a,b,N):                         #funktioniert
    if b >= a:
        h = (b-a)/N
        x = a
        area = (1/3)*(f(a)+f(b))
        i = 1
        while i < N:
            if (i % 2 == 0):
                area += (2/3)*f(x+i*h)
            else:
                area += (4/3)*f(x+i*h)
            i += 1
        area = h*area
    else:
        h = (a-b)/N
        x = b
        area = (1/3)*(f(a)+f(b))
        i = 1
        while i < N:
            if (i % 2 == 0):
                area += (2/3)*f(x+i*h)
            else:
                area += (4/3)*f(x+i*h)
            i += 1
        area = h*area
    return area

#######################################################################################################################################################################################

#Romberg integration ohne Extrapolation, am präsizesten, aber auch am aufwendigsten

def integral_romberg(f,a,b,eps):

    g = 0
    if b < a:
        z = a
        a = b
        b = z
        g = 1
    
    h = b - a 

    x = np.asfarray([h**2])
    y = np.asfarray([])

    Inew = 0
                                        
    I0 = (h/2)*(f(a)+f(b))
    y = np.append(y,I0)

    In = I0
    In1 = 0.5*(I0 + h*f(a+0.5*b))

    i = 2
    o = 0
    z = 0
    
    while o < 1:

        if np.abs((In1-In)) > np.abs(eps*In):
            k = 0
            sum = 0
            deltan = (h/(2**(i-1)))

            In = In1
            
            while k <= (2**(i-1)-1):
                sum += f(a+(k+0.5)*deltan)
                k += 1

            Inew = deltan*sum

            In1 = 0.5*(In+Inew)

            y = np.append(y,In1)
            i += 1
        else:
            o += 1
            
    for j in range(i-2):
        x = np.append(x,(h/(2*(j+1)))**2)

    I_final = y[np.size(y)-1]

    if g == 1:
        return -I_final,y

    return I_final,y
    
#######################################################################################################################################################################################

#Romberg integration mit Neville-Extrapolation (hier durch order gegeben. order 1 entspricht linear usw., order 1 ist am präsizesten), effizienter als ohne

def integral_romberg_extrapol(f,a,b,eps,order):

    g = 0
    if b < a:
        z = a
        a = b
        b = z
        g = 1

    h = b - a 

    x = np.asfarray([h**2])
    y = np.asfarray([])

    Inew = 0
                                        
    I0 = (h/2)*(f(a)+f(b))
    y = np.append(y,I0)

    In = I0
    In1 = 0.5*(I0 + h*f(a+0.5*b))

    i = 2
    

    while i < 3+order:
        k = 0
        sum = 0
        deltan = (h/(2**(i-1)))

        In = In1
        
        while k <= (2**(i-1)-1):
            sum += f(a+(k+0.5)*deltan)
            k += 1

        Inew = deltan*sum

        In1 = 0.5*(In+Inew)

        y = np.append(y,In1)
        x = np.append(x,(h/(2**i))**2)
        i += 1


    x1 = np.delete(x,0)
    y1 = np.delete(y,0)

    x = np.delete(x, np.size(x)-1)
    y = np.delete(y, np.size(y)-1)

    o = 0
    while o < 1:

        if np.abs(neville_algorithm(0,order,x1,y1)-neville_algorithm(0,order,x,y)) > np.abs(eps*neville_algorithm(0,order,x1,y1)):
            k = 0
            sum = 0
            deltan = (h/(2**(i-1)))

            In = In1
            
            while k <= (2**(i-1)-1):
                sum += f(a+(k+0.5)*deltan)
                k += 1

            Inew = deltan*sum

            In1 = 0.5*(In+Inew)

            x = x1
            y = y1
            y1 = np.append(y1,In1)
            x1 = np.append(x1,(h/(2**i))**2)
            x1 = np.delete(x1,0)
            y1 = np.delete(y1,0)

            i += 1
        else:
            if g == 1:
                return -y1[np.size(y1)-1],i
            else:
                return y1[np.size(y1)-1],i

#######################################################################################################################################################################################

def secante(x0,x1,f,precision):
    counter = 0

    a0 = x0
    a1 = x1

    i  = 0

    while i < 1:
        
        if np.abs(f(a1)) > precision:
            a0  = ((f(a1))/(f(a0)-f(a1)))*(a1-a0) + a1
            counter += 1
            if np.abs(f(a0)) > precision:
                a1  = ((f(a0))/(f(a1)-f(a0)))*(a0-a1) + a0
                counter += 1  
            else:
                return counter, a0
        else:
            return counter, a1

#######################################################################################################################################################################################

def regula_falsi(x0,x1,f,precision):
    counter = 0

    a = x0
    b = x1

    anew = 0

    i  = 0

    while i < 1:
        
        if np.abs(f(anew)) > precision:
            anew  = ((f(b))/(f(a)-f(b)))*(b-a) + b
            counter += 1
            if f(b)*f(anew) < 0:
                b = anew
                a = a
            else:
                b = b
                a = anew
        else:
            return counter, anew
        
#######################################################################################################################################################################################     

def bisection(x0,x1,f,precision):
    counter = 0
    eps = np.abs(x1 - x0)
    i = 0
    a = x0
    b = x1
    while i < 1:
        if np.abs(f(a)) > precision and np.abs(f(b)) > precision:
            if f(a)*f(a+0.5*eps) < 0:
                a = a
                b = a + 0.5*eps
                eps = eps/2
                counter += 1
            else:
                a = a + 0.5*eps
                b = b
                eps = eps/2
                counter += 1
        else:
            if np.abs(f(a)) > precision:
                return counter, a
            else:
                return counter, b

#######################################################################################################################################################################################   
    
def QL(A,step):

    dim = np.shape(A)[0]

    # sum = 0                           just to check, whether A stays symmetric, yes it does

    for o in range(step):

        Q = np.identity(dim)

        for i in range(dim-1):
            q = dim - 1 - i
            p = q - 1

            t = -((A[p][q])/(A[q][q]))
            #c = A[q][q]/((A[p][q]**2+A[q][q]**2)**0.5)
            #s = -A[q][p]/((A[q][p]**2+A[q][q]**2)**0.5)
            c = 1/(((t**2)+1)**0.5)
            s = t*c
            P = np.identity(dim)
            P[p][p] = c
            P[q][q] = c
            P[p][q] = s
            P[q][p] = -s

            if i == 0:
                Q = transpose(P)
                continue

            Pt = transpose(P)
            Q = matrix_dot(Q, Pt)

        Qt = transpose(Q)

        L = matrix_dot(Qt,A)

        A = matrix_dot(L,Q)

    return A

#######################################################################################################################################################################################   

def transpose(A):
    row = np.shape(A)[0]
    col = np.shape(A)[1]
    row, col = col, row
    At = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            At[i][j] = A[j][i]
    return At

#######################################################################################################################################################################################
