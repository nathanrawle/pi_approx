"""
Python3 program to calculate Pi using python long integers, binary
splitting and the Chudnovsky algorithm

See: http://www.craig-wood.com/nick/articles/pi-chudnovsky/ for more
info

Nick Craig-Wood <nick@craig-wood.com>
"""

import math
from gmpy2 import mpz, isqrt
from time import time

def pi_chudnovsky_bs(digits):
    """
    Compute int(pi * 10**digits)

    This is done using Chudnovsky's series with binary splitting
    """
    C = 640320
    C3_OVER_24 = C**3 // 24
    def bs(a, b):
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        if b - a == 1:
            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                Pab = Qab = mpz(1)
            else:
                Pab = mpz((6*a-5)*(2*a-1)*(6*a-1))
                Qab = mpz(a*a*a*C3_OVER_24)
            Tab = Pab * (13591409 + 545140134*a) # a(a) * p(a)
            if a & 1:
                Tab = -Tab
        else:
            # Recursively compute P(a,b), Q(a,b) and T(a,b)
            # m is the midpoint of a and b
            m = (a + b) // 2
            # Recursively calculate P(a,m), Q(a,m) and T(a,m)
            Pam, Qam, Tam = bs(a, m)
            # Recursively calculate P(m,b), Q(m,b) and T(m,b)
            Pmb, Qmb, Tmb = bs(m, b)
            # Now combine
            Pab = Pam * Pmb
            Qab = Qam * Qmb
            Tab = Qmb * Tam + Pam * Tmb
        return Pab, Qab, Tab
    # how many terms to compute
    DIGITS_PER_TERM = math.log10(C3_OVER_24/6/2/6)
    N = int(digits/DIGITS_PER_TERM + 1)
    # Calclate P(0,N) and Q(0,N)
    P, Q, T = bs(0, N)
    one_squared = mpz(10)**(2*digits)
    sqrtC = isqrt(10005*one_squared)
    return (Q*426880*sqrtC) // T

if __name__ == '__main__':
    pi_1m = pi_chudnovsky_bs(10**6)
    with open('pi_1m.txt', 'w') as f:
        f.write(str(pi_1m))
    del pi_1m

    pi_100m = pi_chudnovsky_bs(10**8)
    with open('pi_100m.txt', 'w') as f:
        f.write(str(pi_100m))
    del pi_100m