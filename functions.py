class Vector:
    def __init__(self,name):
        self.dummy = []
        self.actualcoord = []
        for varnum in range(nvar):
            self.dummy.append(symbols(name + '_' + str(varnum)))
            self.actualcoord.append(symbols(name + '_' + str(varnum)))
        self.dummy = Matrix(self.dummy)
            
class matrix:
    def __init__(self,name, ncol=nvar, refmatrix=zeros(nvar)):
        self.dummy = zeros(ncol)
        if refmatrix!=zeros(nvar):
            self.actualcoord = refmatrix
        for row in range(ncol):
            for col in range(ncol):
                if refmatrix[row,col]!=0:
                    self.dummy[row,col] = symbols(name + '_' + str(row+1) + '^' + str(col+1))

def isfloat(value):
  try:
    float(eval(value))
    return True
  except:
    return False

def i_n(l,arg,derivative=False):
    z=symbols('z')
    gforward = []
    gforward.append(Pow(z,-1))
    gforward.append(-Pow(z,-2))
    gbackward = []
    auxiliarybackward=-Pow(z,-2)
    gbackward.append(Pow(z,-1))
    for counter in range(1,l+2):
        if counter==1:
            gbackward.append(Add(auxiliarybackward,Mul(-2*(counter-1)+1,Pow(z,-1),gbackward[counter-1])))
        else:
            gbackward.append(Add(gbackward[counter-2],Mul(-2*(counter-1)+1,Pow(z,-1),gbackward[counter-1])))
    if l>1:
        for counter in range(2,l+1):
            gforward.append(Add(gforward[counter-2],Mul(-1,2*(counter-1)+1,Pow(z,-1),gforward[counter-1])))
    modified_spherical_bessel=Add(Mul(gforward[l],sinh(z)),Mul(gbackward[l+1],cosh(z)))
    if derivative==False:
        return simplify(modified_spherical_bessel.subs(z,arg))
    else:
        return simplify(diff(modified_spherical_bessel,z).subs(z,arg))

def legendre(l, m, arg):
    x = symbols('x')
    leg = Pow(Add(Pow(x,2), -1), l)
    for dernum in range(l+m):
        leg = diff(leg, x)
    leg = Mul(Pow(-1, m), Pow(2, -l), Pow(math.factorial(l), -1),
              Pow(Add(1, Mul(-1, Pow(x, 2))), Mul(m, Pow(2, -1))), leg)
    return simplify(leg.subs(x,arg))

def crossorderapplied(u,parameter,equilibrium):
    SS = zeros(nvar,1)
    firstordereval = firstorderderivatives
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            firstordereval[counter1] = Matrix(firstordereval[counter1]).subs(var[nvar+counter2], equilibrium[counter2])
    for counter1 in range(nvar):
        SS = Add(SS, Mul(u[counter1], diff(firstordereval[counter1], parameter)))
    return SS
    
def secondorderapplied(u,v):
    DS = zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            DS = Add(DS, Mul(u[counter1], v[counter2],
                             secondorderderivatives[counter1][counter2]))
    return Mul(Pow(2, -1), DS)
    
def thirdorderapplied(u,v,w):
    TS=zeros(nvar,1)
    for counter1 in range(nvar):
        for counter2 in range(nvar):
            for counter3 in range(nvar):
                TS = Add(TS, Mul(u[counter1], v[counter2], w[counter3],
                                 thirdorderderivatives[counter1][counter2][counter3]))
    return Mul(Pow(6, -1), TS)

def dummyvareval(vector,negativeRHS,coefmat):
    for row in range(nvar):
        vector = vector.subs(negativeRHS.dummy[row], negativeRHS.actualcoord[row])
        for col in range(nvar):
            vector = vector.subs(coefmat.dummy[row,col], coefmat.actualcoord[row,col])
    return vector

def linearsolver(vector,negativeRHS,coefmat):
    vector.actualcoord = linsolve(Add(Mul(coefmat.dummy, vector.dummy), negativeRHS.dummy), list(vector.dummy))
    vector.actualcoord = transpose(Matrix(list(vector.actualcoord)))
    vector.actualcoord = dummyvareval(vector.actualcoord, negativeRHS, coefmat)
    return vector.actualcoord

def add_2_up_to_m(l, m):
    set_check = []
    non_repeated_list = []
    for q1 in range(-l, l+1):
        q2 = m - q1
        if -l <= q2 <=l and {q1, q2} not in set_check:
            set_check.append({q1, q2})
            non_repeated_list.append([q1, q2])
    return non_repeated_list

def add_3_up_to_m(l, m):
    set_check = []
    non_repeated_list = []
    for q1 in range(-l, l+1):
        for q2 in range(-l, l+1):
            q3 = m - q1 - q2
            if -l <= q3 <=l and {q1, q2, q3} not in set_check:
                set_check.append({q1, q2, q3})
                non_repeated_list.append([q1, q2, q3])
    return non_repeated_list

def I_l():
    integral = 0
    for varnum in range(nvar):
        integral = Add(integral, Mul(Pow(Sl[varnum, varnum], 2), phiNF.dummy[varnum], psiNF.dummy[varnum],
                                     Pow(omega[varnum, varnum], -3),
                                     integrate(Mul(Pow(ilomegar.actualcoord[varnum, varnum].
                                                       subs(r, Mul(s, Pow(omega[varnum, varnum], -1))), 2),
                                                   Pow(s,2)), (s, 0, Mul(omega[varnum, varnum], R)))))
        
    integral = Add(integral, Mul(Pow(R,2), phiNF.dummy.dot(psiNF.dummy)))
    
    return integral

def d_p(lNF, mNF, p, q_1, q_2):
    dp = Mul(2*p + 1, math.factorial(p - abs(q_1 + q_2)),
             Pow(2, -1), Pow(math.factorial(p + abs(q_1 + q_2)), -1),
             integrate(Mul(legendre(lNF, abs(q_1), x),
                           legendre(lNF, abs(q_2), x),
                           legendre(p, abs(q_1 + q_2), x)), (x, -1, 1)))
    return dp

def first_term_C3(lNF, mNF, q_1, q_2, q_3):
    summation = 0
    for p in range(abs(q_1 + q_2), 2*lNF + 1):
        summation = Add(summation, Mul(DS_phiup2[p - abs(q_1 + q_2)],
                                       integrate(Mul(legendre(lNF, abs(q_3), x),
                                                     legendre(p, abs(q_1 + q_2), x),
                                                     legendre(lNF, abs(mNF), x)), (x, -1, 1))))
    
    return summation
        
def fourth_integral(lNF, mNF, q_1, q_2, q_3):
    integral = integrate(Mul(legendre(lNF, abs(q_1), x),
                  legendre(lNF, abs(q_2), x),
                  legendre(lNF, abs(q_3), x),
                  legendre(lNF, abs(mNF), x)), (x, -1, 1))
    
    return integral

def non_repeated_element(lista):
    if len({lista[0], lista[1]})==1:
        return 2
    elif len({lista[0], lista[2]})==1:
        return 1
    else:
        return 0