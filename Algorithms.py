from math import log2

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

def fast_convolution(f, g, field):
    

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