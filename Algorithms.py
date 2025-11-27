from math import log2
from copy import deepcopy

def karatsuba(p1, p2, field):
    #print("Karatsuba")
    #print_poly(p1)
    #print_poly(p2)
    maxdeg = max(len(p1), len(p2)) - 1
    r = log2(maxdeg + 1)
    if r != int(r):
        r = int(r) + 1
    r = int(r)
    if r == 0:
        return [(p1[0] * p2[0]) % field]
    while len(p1) < 2**r:
        p1 = p1 + [0]
    while len(p2) < 2**r:
        p2 = p2 + [0]
    f0 = p1[:2**(r-1)]
    f1 = p1[2**(r-1):]
    g0 = p2[:2**(r-1)]
    g1 = p2[2**(r-1):]
    #print_poly(f0)
    #print_poly(f1)
    #print_poly(g0)
    #print_poly(g1)
    f0_plus_f1 = poly_add(f0, f1, field)
    g0_plus_g1 = poly_add(g0, g1, field)
    while len(f0_plus_f1) < 2**(r-1):
        f0_plus_f1 = f0_plus_f1 + [0]
    while len(g0_plus_g1) < 2**(r-1):
        g0_plus_g1 = g0_plus_g1 + [0]
    f0g0 = karatsuba(f0, g0, field)
    f1g1 = karatsuba(f1, g1, field)
    f0pf1g0pg1 = karatsuba(f0_plus_f1, g0_plus_g1, field)
    star2 = [0]*(2**(r-1)) + poly_minus(poly_minus(f0pf1g0pg1, f0g0, field), f1g1, field)
    result = poly_add(poly_add(f0g0, [0]*(2**r) + f1g1, field), star2, field)
    result = enforce_poly(result, field)
    return result


def print_poly(p):
    string = ""
    for i in range(len(p)):
        if p[i] == 0:
            continue
        x_string = "x^" + str(i)
        if i == 0:
            x_string = ""
        coeff_string = str(p[i])
        if p[i] == 1 and i != 0:
            coeff_string = ""
        string = string + coeff_string + x_string
        if i < len(p) - 1:
            string = string + " + "
    print(string)

def enforce_poly(p, field):
    if len(p) == 0:
        return p
    while p[-1] == 0:
        p = p[:-1]
        if len(p) == 0:
            return p
    for i in range(len(p)):
        p[i] = p[i] % field
    return p


def poly_mult(p1, p2, field):
    p3 = [0 for i in range(len(p1) + len(p2) - 1)]
    for i in range(len(p1)):
        for j in range(len(p2)):
            coeff = (p1[i] * p2[j]) % field
            p3[i+j] = (p3[i+j] + coeff) % field
    p3 = enforce_poly(p3, field)
    return p3

def poly_add(p1, p2, field):
    if len(p1) < len(p2):
        p1 = p1 + [0]*(len(p2) - len(p1))
    elif len(p2) < len(p1):
        p2 = p2 + [0]*(len(p1) - len(p2))
    p3 = [(p1[i] + p2[i])%field for i in range(max(len(p1), len(p2)))]
    p1 = enforce_poly(p1, field)
    p2 = enforce_poly(p2, field)
    p3 = enforce_poly(p3, field)
    return p3

def poly_minus(p1, p2, field):
    
    if len(p1) < len(p2):
        p1 = p1 + [0]*(len(p2) - len(p1))
    elif len(p2) < len(p1):
        p2 = p2 + [0]*(len(p1) - len(p2))
    p3 = [(p1[i] - p2[i])%field for i in range(max(len(p1), len(p2)))]
    p1 = enforce_poly(p1, field)
    p2 = enforce_poly(p2, field)
    p3 = enforce_poly(p3, field)
    return p3

def inverse(x, field):
    for i in range(field):
        if (i * x) % field == 1:
            return i
        
def NewtonIteration(poly, l, field):
    g_prev = [1]
    r = log2(l)
    if r != int(r):
        r = int(r) + 1
    for i in range(1,r+1):
        g = poly_minus(poly_mult(g_prev, [2], field), poly_mult(poly, poly_mult(g_prev, g_prev, field), field), field)
        if len(g) <= 2**i:
            g_prev = g
            continue
        g_prev = g[:2**i]
    return enforce_poly(g_prev[:l], field)

def fast_fourier_transform(f, omegas, field):
    maxdeg = len(f) - 1
    r = log2(maxdeg + 1)
    if r != int(r):
        r = int(r) + 1
    r = int(r)
    while len(f) < 2**r:
        f = f + [0]
    if r == 0:
        return [f[0]]
    f0 = f[:2**(r-1)]
    f1 = f[2**(r-1):]
    r0 = poly_add(f0, f1, field)
    omega = omegas[1]
    r1_semi = poly_minus(f0, f1, field)
    r1 = [((omega**i)*r1_semi[i])%field for i in range(len(r1_semi))]
    fft_r0 = fast_fourier_transform(r0, omegas[::2], field)
    fft_r1 = fast_fourier_transform(r1, omegas[::2], field)
    result = [x for i in range(len(fft_r0)) for x in (fft_r0[i], fft_r1[i])]
    return result

def eval_poly(f, x, field):
    return sum([((x**i)*f[i])%field for i in range(len(f))])%field

def fast_convolution(f, g, omega, n, field):
    omegas = [(omega**i)%field for i in range(n)]
    alpha = fast_fourier_transform(f, omegas, field)
    beta = fast_fourier_transform(g, omegas, field)
    gamma = [(alpha[i]*beta[i])%field for i in range(len(alpha))]
    omega_inv = inverse(omega, field)
    omegas_inv = [(omega_inv**i)%field for i in range(n)]
    dft_res = fast_fourier_transform(gamma, omega_inv, field)
    n_inv = inverse(n, field)
    result = [(dft_res[i]*n_inv)%field for i in range(len(dft_res))]
    return result

class MultiVarPoly:
    def __init__(self, degrees, coefficients, field, numvars):
        self.degrees = degrees
        self.coefficients = coefficients
        self.field = field
        self.vars = numvars

    def __str__(self):
        string = ""
        for i in range(len(self.degrees)):
            degi = self.degrees[i]
            x_string = ""
            for j in range(len(degi)):
                if degi[j] > 0:
                    x_string = x_string + "x_{" + str(j+1) + "}"
                    if degi[j] > 1:
                        x_string = x_string + "^" + str(degi[j])
            if self.coefficients[i] != 1 or x_string == "":
                string = string + str(self.coefficients[i]) + x_string
            else:
                string = string + x_string
            if i < len(self.degrees) - 1 and (self.coefficients[i+1] > 0):
                string = string + " + "
            elif i < len(self.degrees):
                string = string + " "
        return string

    def mdeglessthanorequal(self, mdeg1, mdeg2):
        for i in range(len(mdeg1)):
            if mdeg1[i] < mdeg2[i]:
                return 1
            if mdeg2[i] < mdeg1[i]:
                return -1
        return 0

    def enforce_structure(self):
        dict_struct = {}
        for i in range(len(self.degrees)):
            degtuple = tuple(self.degrees[i])
            if degtuple in dict_struct.keys():
                dict_struct[degtuple] += self.coefficients[i]
                if str(self.field) != "Q":
                    dict_struct[degtuple] = dict_struct[degtuple] % self.field
                if abs(dict_struct[degtuple]) < 10**(-14):
                    del dict_struct[degtuple]
            elif abs(self.coefficients[i]) > 10**(-14):
                dict_struct[degtuple] = self.coefficients[i]
        
        listdegs = list(dict_struct)
        listcoeffs = [dict_struct[e] for e in listdegs]
        listdegs2 = [list(e) for e in listdegs]
        #print("postliststuff")
        #print(listdegs2)
        sorteddegs, sortedcoeffs = self.sort(listdegs2, listcoeffs)
        if len(sorteddegs) == 0:
            sorteddegs = [[0 for _ in range(self.vars)]]
            sortedcoeffs = [0]
        self.degrees = sorteddegs
        self.coefficients = sortedcoeffs

    def sort(self, degs, coeffs):
        if len(degs) < 2:
            return degs, coeffs
        split = len(degs)//2
        degs1 = degs[:split]
        coeffs1 = coeffs[:split]
        degs2 = degs[split:]
        coeffs2 = coeffs[split:]
        sorteddegs1, sortedcoeffs1 = self.sort(degs1, coeffs1)
        sorteddegs2, sortedcoeffs2 = self.sort(degs2, coeffs2)
        sorteddegs3, sortedcoeffs3 = self.merge(sorteddegs1, sortedcoeffs1, sorteddegs2, sortedcoeffs2)
        return sorteddegs3, sortedcoeffs3

    def merge(self, degs1, coeffs1, degs2, coeffs2):
        degs3 = []
        coeffs3 = []
        i1 = 0
        i2 = 0
        while i1 < len(degs1) and i2 < len(degs2):
            degi1 = degs1[i1]
            degi2 = degs2[i2]
            leq = self.mdeglessthanorequal(degi1, degi2)
            if leq == 1:
                degs3.append(degi1)
                coeffs3.append(coeffs1[i1])
                i1 += 1
            else:
                degs3.append(degi2)
                coeffs3.append(coeffs2[i2])
                i2 += 1
        while i1 < len(degs1):
            degs3.append(degs1[i1])
            coeffs3.append(coeffs1[i1])
            i1 += 1
        while i2 < len(degs2):
            degs3.append(degs2[i2])
            coeffs3.append(coeffs2[i2])
            i2 += 1
        return degs3, coeffs3
    
    def __add__(self, other):
        p3degrees = []
        p3coeffs = []
        i1 = 0
        i2 = 0
        while i1 < len(self.degrees) and i2 < len(other.degrees):
            mdegi1 = self.degrees[i1]
            mdegi2 = other.degrees[i2]
            leq = self.mdeglessthanorequal(mdegi1, mdegi2)
            if leq == 0:
                coeff = 0
                if str(self.field) != "Q":
                    coeff = (self.coefficients[i1] + other.coefficients[i2]) % field
                else:
                    coeff = (self.coefficients[i1] + other.coefficients[i2])
                if coeff != 0:
                    p3degrees.append(mdegi1)
                    p3coeffs.append(coeff)
                i1 += 1
                i2 += 1
            elif leq == 1:
                p3degrees.append(mdegi1)
                p3coeffs.append(self.coefficients[i1])
                i1 += 1
            elif leq == -1:
                p3degrees.append(mdegi2)
                p3coeffs.append(other.coefficients[i2])
                i2 += 1
        while i1 < len(self.degrees):
            p3degrees.append(self.degrees[i1])
            p3coeffs.append(self.coefficients[i1])
            i1 += 1
        while i2 < len(other.degrees):
            p3degrees.append(other.degrees[i2])
            p3coeffs.append(other.coefficients[i2])
            i2 += 1
        p3 = MultiVarPoly(p3degrees, p3coeffs, self.field, self.vars)
        p3.enforce_structure()
        return p3

    def __mul__(self, other):
        p3degrees = []
        p3coeffs = []
        for i in range(len(self.degrees)):
            for j in range(len(other.degrees)):
                deg1 = self.degrees[i]
                coeff1 = self.coefficients[i]
                deg2 = other.degrees[j]
                coeff2 = other.coefficients[j]
                deg3 = [deg1[k] + deg2[k] for k in range(len(deg1))]
                coeff3 = 0
                if str(self.field) == "Q":
                    coeff3 = coeff1 * coeff2
                else:
                    coeff3 = (coeff1 * coeff2)%self.field
                p3degrees.append(deg3)
                p3coeffs.append(coeff3)
        p3 = MultiVarPoly(p3degrees, p3coeffs, self.field, self.vars)
        #print("postmult")
        #print(p3)
        p3.enforce_structure()
        return p3

    def __sub__(self, other):
        p3 = MultiVarPoly([[0 for _ in range(len(self.degrees[0]))]], [-1], self.field, self.vars) * other
        return self + p3

    def lt(self):
        return MultiVarPoly([self.degrees[-1]], [self.coefficients[-1]], self.field, self.vars)
    
    def lm(self):
        return MultiVarPoly([self.degrees[-1]], [1], self.field, self.vars)

    def lc(self):
        return self.coefficients[-1]

    def mdeg(self):
        return self.degrees[-1]

    def divisible(self, other):
        mdeg1 = self.mdeg()
        mdeg2 = other.mdeg()
        for i in range(len(mdeg1)):
            if mdeg2[i] > mdeg1[i]:
                return False
        return True

def Syzygy(p1, p2):
    lt1 = p1.lt()
    lt2 = p2.lt()
    maxdeg = [max(lt1.degrees[0][i], lt2.degrees[0][i]) for i in range(len(lt1.degrees[0]))]
    coeffs1 = 0
    coeffs2 = 0
    if str(p1.field) == "Q":
        coeffs1 = 1/lt1.lc()
        coeffs2 = 1/lt2.lc()
    else:
        coeffs1 = inverse(lt1.lc(), p1.field)
        coeffs2 = inverse(lt2.lc(), p1.field)
    
    part1 = MultiVarPoly([[maxdeg[i] - lt1.degrees[0][i] for i in range(len(maxdeg))]], [coeffs1], p1.field, p1.vars)
    part2 = MultiVarPoly([[maxdeg[i] - lt2.degrees[0][i] for i in range(len(maxdeg))]], [coeffs2], p1.field, p1.vars)
    #print("hej")
    #print(part1 * p1)
    #print(part2 * p2)
    return (part1 * p1) - (part2 * p2)

def multiVarRemainder(p1, G):
    r = MultiVarPoly([[0 for i in range(p1.vars)]], [0], p1.field, p1.vars)
    Q = [MultiVarPoly([[0 for i in range(p1.vars)]], [0], p1.field, p1.vars) for _ in range(len(G))]
    f = deepcopy(p1)
    while len(f.coefficients) != 1 or f.coefficients[0] != 0:
        divisor = None
        divindex = 0
        for i, g in enumerate(G):
            if f.divisible(g):
                divisor = g
                divindex = i
                break
        if divisor == None:
            r = r + f.lt()
            f = f - f.lt()
            r.enforce_structure()
            f.enforce_structure()
        else:
            coeff = 0
            if str(p1.field) != "Q":
                coeff = (f.lc() * inverse(divisor.lc(), p1.field)) % p1.field
            else:
                coeff = f.lc() / divisor.lc()
            quotient = MultiVarPoly([[f.mdeg()[i] - divisor.mdeg()[i] for i in range(len(f.mdeg()))]], [coeff], p1.field, p1.vars)
            quotient.enforce_structure()
            Q[divindex] = Q[divindex] + quotient
            f = f - (quotient * divisor)
    return Q, r

def Buchbergers(G):
    while True:
        S = []
        print("start loop")
        for i in range(len(G) - 1):
            for j in range(i+1, len(G)):
                syz = Syzygy(G[i], G[j])
                print(syz)
                _, r = multiVarRemainder(syz, G)
                if r.coefficients[0] != 0:
                    S.append(r)
        if len(S) == 0:
            return G
        print("iter")
        for e in S:
            print(e)
        G = G + S

if __name__ == "__main__":
    field = 5
    f = [1,0,1]
    l = 5
    finv = NewtonIteration(f, l, field)
    print_poly(finv)
    print_poly(poly_mult(f, finv, field))
    f_finv = karatsuba(f, finv, field)
    print_poly(f_finv)
    print(f_finv)
    print("DFT Testing")
    field = 41
    omega = 14
    omegas = [(omega**i)%field for i in range(8)]
    f = [6, 2, 0, 0, 3, 0, 2, 1]
    dft_brute_force = [eval_poly(f, o, field) for o in omegas]
    print(dft_brute_force)
    print(fast_fourier_transform(f, omegas, field))

    print("Multipoly testing")
    p1 = MultiVarPoly([[0,0,0], [0,0,2], [0,2,0]], [-1, 1, 1], "Q", 3)
    print(p1)
    p2 = MultiVarPoly([[0,0,1], [0,1,1], [1,1,0]], [-2, 2, -2], "Q", 3)
    print(p2)
    p3 = p1 + p2
    print(p3)
    p4 = p1 * p2
    print(p4)
    p5 = Syzygy(p1, p2)
    print(p5)
    print(p1 - p1)
    print("Remainder testing")
    g1 = MultiVarPoly([[0,0,0], [0,0,2], [0,2,0]], [-1, 1, 1], "Q", 3)
    g2 = MultiVarPoly([[0,0,1], [0,1,1], [1,1,0]], [-2, 2, -2], "Q", 3)
    g3 = MultiVarPoly([[0,0,0], [0,1,0], [0,2,0], [1,0,1]], [1, -2, 1, -2], "Q", 3)
    G = [g1, g2, g3]
    f = Syzygy(g1, g2)
    Qs, r = multiVarRemainder(f, G)
    for e in Qs:
        print(e)
    print(r)
    print("grobner testing")
    grobner = Buchbergers(G)
    print("basis1")
    for e in grobner:
        print(e)
        pass
    grobner = grobner[:-7]
    for e in grobner:
        #print(e)
        pass
    grobner = [grobner[3], grobner[6], grobner[8]]
    for e in grobner:
        print(e)
    mingrobner = [MultiVarPoly([[0,0,0]], [-1], "Q", 3)*grobner[0], MultiVarPoly([[0,0,0]], [1/(grobner[1].lc())], "Q", 3)*grobner[1], MultiVarPoly([[0,0,0]], [1/(grobner[2].lc())], "Q", 3)*grobner[2]]
    for e in mingrobner:
        print(e)