from sympy import *
from sympy.abc import x, alpha
from sympy import lcm
import numpy as np
import random
import copy
from itertools import *
import galois


def order(a,p):
    cur = a%p
    counter =0
    while True:
        counter +=1
        if (cur==1):
            return counter
        cur = (cur*a)%p
    raise Exception("Couldn't determine order of poly")

def fastExp(poly, exp, irr_poly):
    if exp == 0:
        return Poly(1, x, domain = poly.domain)
    if exp == 1:
        return poly%irr_poly
    res = fastExp(poly, exp//2, irr_poly)
    res = (res*res)%irr_poly
    if exp%2==1:
        res = (res * poly)%irr_poly
    return res
    
def findNthRoot(q, m, irr_poly, n):
    assert((q**m -1)%n==0)
    for i in range(2*(q**m)):
        coeffs = [0] + [random.randint(0,q-1) for _ in range(m)]
        if sum(coeffs)==0:
            continue
        poly = Poly(coeffs, x, domain = GF(q))
        r = order(poly, irr_poly)
        if r%n==0:
            return fastExp(poly, (r//n), irr_poly)
    raise Exception("Couldn't determine nth root")

def evaluatePoly(poly, min_poly, q, m, irr_poly):#poly has variable x, min_poly has alpha
    res = Poly(0, x, domain = GF(q))
    for i in min_poly.all_terms():
        if i[1]==0:
            continue
        res = (res + i[1]*fastExp(poly, i[0][0], irr_poly))%irr_poly
    return res

# return minimum polynomial of poly in GF(q^m) with reducing polynomial irr_poly
#may have negative coeffs
def minPoly(poly, q, m, irr_poly):#poly and irr_poly are both in variable x
    for i in range(1,m+1):
        if m%i!=0:
            continue
        for coeffs in product(range(q), repeat = i):
            candidate_poly = Poly([1]+list(coeffs), alpha, domain = GF(q))
            if evaluatePoly(poly, candidate_poly, q, m, irr_poly).is_zero:
                return candidate_poly
    raise Exception("Couldn't find minimal polynomial")

def polyToInt(poly, q):
    return sum(int((c%q)*(q**i)) for i,c in enumerate(reversed(poly.all_coeffs())))

def intToPoly(k, q):
    deg = 0
    res = Poly(0,x,domain = GF(q))
    while(k>0):
        res = res + Poly((k%q)*x**deg, x, domain = GF(q))
        deg = deg+1
        k = k//q
    return res


class BCHCode:
    def __init__(self, q, n, d, c):
        assert(isprime(q))
        assert (np.gcd(q, n)==1)
        assert(d<n), "The value of d is too large."
        self.q = q
        self.n = n
        self.d = d
        self.c = c
        self.m = order(q,n)
        self.t = (d-1)//2
        self.elem = []
        self.gf = galois.GF(self.q,self.m)
        (self.g_poly,self.irr_poly, self.nth_root) = self.__gen()
    
    def __gen(self):
        coeffs = [int(c) for c in self.gf.irreducible_poly.coeffs]
        irr_poly = Poly(coeffs, x, domain = GF(self.q))
        nth_root = findNthRoot(self.q, self.m, irr_poly, self.n)
        self.__calcElem(nth_root, irr_poly)
        g_poly = Poly(1, alpha, domain = GF(self.q))
        current_root =  fastExp(nth_root, self.c, irr_poly)
        for i in range(self.c, self.c + self.d -1):
            g_poly = lcm(g_poly, minPoly(current_root, self.q, self.m, irr_poly))
            current_root = (current_root*nth_root)%irr_poly
        return(g_poly, irr_poly, nth_root)

    def __calcElem(self, nth_root, irr_poly):#called once during _gen
        cur = Poly(1, x, domain = GF(self.q))
        for i in range(self.n):
            self.elem.append(cur)
            cur = (cur * nth_root)%irr_poly
        
    #message is a q-ary list of any size at most n - g_poly.degree. Automatically padded with zeros to the right
    def encode(self, message):
        assert(len(message)+self.g_poly.degree(alpha)<=self.n), "Message too long!"
        m_poly = Poly(message, alpha, domain = GF(self.q))*Poly(alpha**(self.n - len(message)), alpha, domain = GF(self.q))
        r_poly = m_poly % self.g_poly
        m_poly = m_poly-r_poly
        code = [0]*(self.n - m_poly.degree(alpha)-1)
        for i in m_poly.all_coeffs():
            code.append(i%self.q)
        return code

    #code is a q-ary list of size n
    def decode(self, code):
        code1 = copy.deepcopy(code)
        code_poly = Poly(code1, alpha, domain = GF(self.q)) 
        syndromes = []
        clean = True
        for i in range(self.c, self.c + self.d - 1):
            s = evaluatePoly( self.elem[i%self.n], code_poly, self.q, self.m, self.irr_poly)
            syndromes.append(s)#poly in x variable
            if not s.is_zero:
                clean = False
        if clean:
            print("No errors detected!")
            return code1[:self.n - self.g_poly.degree(alpha)]
        
        S = self.gf([[polyToInt(syndromes[i+j], self.q) for j in range(self.t)] for i in range(self.t)])
        while np.linalg.det(S)==0 and S.size>1:
            S = S[:-1, :-1]

        v = S.shape[0]
        C = self.gf([[polyToInt(syndromes[i+v], self.q)] for i in range(v)])
        print(f"Not clean! Detected {v} errors.")
        Inv = np.linalg.inv(S)
        lambdas = Inv @ C

        error_locator_poly = Poly(-1, x, domain = GF(self.q))#-1 since I solved for -lambda instead
        for i in range(len(lambdas)):           
            error_locator_poly = error_locator_poly + intToPoly(int(lambdas[i][0]),self.q)*(Poly(alpha**(v-i), alpha, domain = GF(self.q)))

        #test each power of nth_root to find the solutions
        error_pos = []
        for j in range(self.n):#try nth_root**j
            res = Poly(0, x, domain = GF(self.q))
            for i in error_locator_poly.terms():
                if i[1]==0:
                    continue
                res = (res + i[1]*fastExp(self.elem[j], i[0][1], self.irr_poly)*Poly(x**i[0][0], x, domain = GF(self.q)))%self.irr_poly
            if res.is_zero:
                print("Error found at position: ", (j-1)%self.n)
                if self.q ==2:
                    code1[(j-1)%self.n] = 1 if code1[(j-1)%self.n] ==0 else 0
                else:
                    error_pos.append((self.n-j)%self.n)
        if self.q==2:
            return code1[:self.n - self.g_poly.degree(alpha)]
        
        S = self.gf([[polyToInt( fastExp(self.elem[error_pos[j]], self.c+i, self.irr_poly), self.q) for j in range(v)] for i in range(v)])
        C = self.gf([[polyToInt(syndromes[i], self.q)] for i in range(v)])
        Inv = np.linalg.inv(S)
        errors = Inv @ C

        for j in range(v):
            code1[(self.n-error_pos[j]-1)%self.n] = (code1[(self.n-error_pos[j]-1)%self.n]-int(errors[j][0]))%self.q
        return code1[:self.n - self.g_poly.degree(alpha)]
