#!/usr/bin/env python
# coding: utf-8

# # Different code implementations of various algorithms for approximating Pi
# Following along with the walk-through by Nick Craig-Wood, starting [here](https://www.craig-wood.com/nick/articles/pi-gregorys-series/)

# In[1]:


import math


# In[2]:


from decimal import *


# In[3]:


D = Decimal


# In[4]:


from time import time


# ## Gregory's Series
# https://www.craig-wood.com/nick/articles/pi-gregorys-series/

# In[5]:


def gregory_pi(n):
    '''
    Implementation of:
    
    pi/4 = 1 - 1/3 + 1/5 - 1/7 + … - 1/2(n-1)+1 + 1/2n+1
    
    aka Gregory's series or the Leibniz formula.
    
    Strength: simple
    Weakness: Not very efficient. 10^8 iterations for only 7 dp precision
    
    Using Decimal() eliminates(?) error introduced by limited precision of float() numbers
    '''
    
    sigma = D(0)
    num = D(1)
    den = D(1)
    
    for i in range(n):
        sigma += num/den
        num = -num
        den = den + 2
    
    return 4*sigma
        


# In[6]:


g_start = time()
g_approx = gregory_pi(100000000)
g_elapsed = time() - g_start

print(g_approx, f'{g_elapsed:.2f}s')


# ## Archimedes method
# https://www.craig-wood.com/nick/articles/pi-archimedes/
# 
# Starting from a square inscribed within a circle with radius 1 (thus circumference 2\*Pi), use Pythagoras to iteratively calculate the edge-length of successive bisections. After n iterations, multiply the last length calculated by the half the number of edges in the last inscribed polygon (2^n + 2) to obtain the estimated half-circumference and thus the estimate of Pi. (This makes more sense with the explanation given at the link above)

# In[7]:


def archimedes_pi(prec):
    '''
    approximate pi as the perimeter of an n-sided polygon inscribed within a circle.
    
    For r=1, as n tends to infinity, the perimeter tends to 2*pi.
    
    Strengths: easy to understand how it works and way faster than Gregory's series
    Weakness: Struggles beyond 1 000dp and slow (square root in loop is expensive)
    '''
    dp = prec+1 #<-- for 'prec' entered to match real d.p. in final result
    old_result = None
    
    # Simple if inefficient way to walk up to the precision desired.
    # Assumes that the max number of iterations needed will be < 10*precision desired
    for n in range(int(1.67*dp), 10*dp):
        # set 'double' precision for calculations
        getcontext().prec = 10*dp
        
        # seed with a square (dn = sqrt(2))
        d_n_squared = D(2)
        half_polygon = 2
        

        for i in range(n):
            d_n_squared = 2 - 2*(1 - d_n_squared/4).sqrt()
            half_polygon *= 2
            
        result = half_polygon * d_n_squared.sqrt()
        
        # set single-precision for result
        getcontext().prec = dp
        result = +result
        if result == old_result:
            return result
            break
        old_result = result


# In[8]:


a_start = time()
a_approx = archimedes_pi(2000)
a_elapsed = time() - a_start
print('Archimedes: ', a_approx, f'Time: {a_elapsed:.2f}s')


# So about 10 minutes for 2000 digits!
# 
# Double-precision float maths works for a few hundred dp, but for 1 000dp I needed to go to triple-precision.
# I couldn't quite get 2 000dp even at 10\*prec (last three decimals should be 009 not 010). Not too bad!
# 
# (validated against https://www.piday.org/million/)

# ## Machin's Formula
# https://www.craig-wood.com/nick/articles/pi-machin/
# 
# Machin's formula uses arctan. Python's built-in math.atan() is too slow. First we define two alternative functions for calculating arctan using fixed-point arithmetic: a standard version and Euler's accelerated version. Then we define Machin's original formula. There are several variations of Machin's formula, the best one being Gauss's, which is also defined below.

# In[9]:


def arctan_1_over(x, one=10**100):
    '''
    custom function for calculating arctan(1/x) with fixed-point arithmetic
    
    arctan(1/x) = 1/x - 1/(3x^3) + 1/(5x^5) - 1/(7x^7) + 1/(9x^9) - … (x >= 1)
    
    '''
    # Seeds
    power = one // x
    total = power
    x_squared = x*x
    divisor = 1
    
    while 1:
        power = -power // x_squared
        divisor += 2
        power += divisor // 2 # round the division (reduces the number of erroneous least-significant digits)
        delta = power // divisor
        if delta == 0:
            break
        total += delta
    return total


# In[10]:


def euler_atan_1_over(x, one=10**100):
    '''
    custom function for calculating arctan(1/x) with fixed-point arithmetic
    euler's accelerated formula for approximating arctan 1/x
    
    arctan(1/x) = x/(1+x^2) + (2/3)*(x/(1+x^2)^2 + (2/3)*(4/5)*(x/(1+x^2)^3) + …
    '''
    
    # seeds
    x_squared = x*x
    x_squared_plus_1 = x_squared + 1
    term = (x * one) // x_squared_plus_1
    total = term
    two_n = 2
    
    while True:
        divisor = (two_n + 1) * x_squared_plus_1
        term *= two_n
        term += divisor // 2
        term = term // divisor
        if term == 0:
            break
        total += term
        two_n += 2
    return total


# In[11]:


def str_format_fix(result):
    '''
    Fixed-point arithmetic returns the digits of pi encoded as an integer with log10(one) digits
    Decimal can't handle precision beyond 999 999 999 999 999 999 (10**18 - 1) on a 64-bit system
    To get something that looks like a decimal, I've cast the output to a string
    '''
    result = str(result)
    return result[0]+'.'+result[1:]


# In[12]:


def machin_pi(arctan, one=10**100):
    '''
    This function calculates pi to log10(one) precision using Machin's Formula
    Machin determined the following fast approximation of pi:
    
    pi/4 = 4*arctan(1/5) - arctan(1/239)
    
    This depends on an arctan function that works using fixed-point arithmetic.
    '''
    
    result = 4*(4*arctan(5, one) - arctan(239, one))
    return str_format_fix(result)


# In[13]:


def gauss_pi(arctan, one=10**100):
    '''
    Gauss's apparently slightly better formula for approximating Pi
    '''
    result = 4*(12*arctan(18, one) + 8*arctan(57, one) - 5*arctan(239, one))
    return str_format_fix(result)


# In[14]:


m_start = time()
m_approx = machin_pi(arctan_1_over, 10**100000)
m_elapsed = time() - m_start

print(f'Machin: {m_approx[0:101]}…{m_approx[-100:]}\nTime: {m_elapsed:.2f}s\n')

m_e_start = time()
m_e_approx = machin_pi(euler_atan_1_over, 10**100000)
m_e_elapsed = time() - m_e_start

print(f'Machin-Euler: {m_e_approx[0:101]}…{m_e_approx[-100:]}\nTime: {m_e_elapsed:.2f}s\n')

g_start = time()
g_approx = gauss_pi(arctan_1_over, 10**100000)
g_elapsed = time() - g_start

print(f'Gauss: {g_approx[0:101]}…{g_approx[-100:]}\nTime: {g_elapsed:.2f}s\n')

g_e_start = time()
g_e_approx = gauss_pi(euler_atan_1_over, 10**100000)
g_e_elapsed = time() - g_e_start

print(f'Gauss-Euler: {g_e_approx[0:101]}…{g_e_approx[-100:]}\nTime: {g_e_elapsed:.2f}s\n')


# For 100 000 digits, using the Euler arctan formula is roughly twice as fast as the basic one on my machine. Gauss makes a modest improvement to the speed of Machin. All can easily calculate 10^5 digits within seconds, and 10^6 in about 30 minutes.
# 
# The last 5 digits should be …24646. I would be interested in knowing whether anything can be done to eliminate the error! For the values of 'one' I've tried, the Machin-type formulas give roughly log10(log10(one)) erroneous least-significant digits.

# ## Chudnovsky Algorithm
# https://www.craig-wood.com/nick/articles/pi-chudnovsky/
# 
# This is the heavyweight of Pi approximation algorithms. Starting from a formula by Ramanujan, brothers David and Gregory Chudnovsky derived this fast algorithm in 1988. As of 2 September 2020, this algorithm has been responsible for the last 8 world-record approximations, including 50Tn digits on 29 Jan 2020 [\[1\]](https://en.wikipedia.org/wiki/Chronology_of_computation_of_%CF%80#With_electronic_computers_(1949%E2%80%93)).
# 
# The built-in math.sqrt method won't handle the size of the numbers generated by the Chudnovsky algorithm, so Newton's method of approximation is used to define a square root function that can be used.

# In[15]:


def sqrt(a, one=10**6):
    '''
    This custom function implements Newton's method applied to refining an initial guess
    at the square-root 'x' of a number 'a' (i.e. x^2 = a), where:
    
    x[n+1] = x[n] - (x^2 - a) / 2x
    
    uses fixed-point arithmetic
    Will break if log(a) << log(one)
    '''
    
    floating_point_precision = 10**16
    a_float = float((a * floating_point_precision)//one) / floating_point_precision
    x = (int(math.sqrt(a_float) * floating_point_precision)*one) // floating_point_precision
    a_one = a * one
    
    while True:
        old_x = x
        x = (x + a_one // x) // 2
        if x == old_x:
            break
    
    return x


# In[16]:


def chudnovsky_pi(one=1000000):
    """
    Calculate pi using the Chudnovskys's series using the fixed-point value for 'one' passed in.
    A computable form of the algorithm is as follows:
    
    Pi = (426880*sqrt(10005)) / (13591409a + 545140134b)
    
    (LET C = 640320)
    
    a = sum_k((-1^k)*(6k)! / ((3k)!)*((k!)^3)*(C^(3k)))
    b = sum_k((-1^k)*(6k)!*k / ((3k)!)*((k!)^3)*(C^(3k)))
    
    a_k = (-1^k)((6k)!) / ((3k)!)((k!)^3)(C^(3k))
    b_k = a_k * k
    a_(k+1) = a_k * (6(k+1)-5)(2(k+1)-1)(6(k+1)-1) / ((k+1)^3)(C^3)
    """
    
    # initialise values of a_k, a_sum, b_sum for k=0 and set-up k for k=1
    k = 1
    a_k = one
    a_sum = one
    b_sum = 0
    
    # store value of (C^3)/24 in memory for efficient arithmetic in the loop
    C = 640320
    C3_OVER_24 = C**3 // 24
    
    # calculate iterative terms of a_k and b_k, and add them to a_sum and b_sum
    while 1:
        a_k *= -(6*k-5)*(2*k-1)*(6*k-1)
        a_k //= k*k*k*C3_OVER_24 #<-- presumably more efficient than k**3
        a_sum += a_k
        b_sum += k * a_k
        k += 1
        # when k^3 >> (A*k-B)^3, integer division gives a_k = 0
        if a_k == 0:
            break
    
    # perform final non-iterative steps to obtain fixed-point value for Pi
    total = 13591409*a_sum + 545140134*b_sum
    pi = (426880*sqrt(10005*one, one)*one) // total
    return str_format_fix(pi)


# In[17]:


c_start = time()
c_approx = chudnovsky_pi(10**10**6)
c_elapsed = time() - c_start

print(f'Chudnovsky: {c_approx[:101]}…{c_approx[-100:]}\nTime: {c_elapsed:.2f}')


# Ten minutes for 1M digits! …almost: last 12 digits should be 105779458151
# 
# Further optimisation is possible using something called 'binary splitting' and the gmpy2 module for precision maths. See the [original walkthrough](https://www.craig-wood.com/nick/articles/pi-chudnovsky/) for details. One or the other seems to eliminate the error in the least significant digits.
# 
# The fully optimised code gives me 1M digits in less than a second, and can provide **100M digits in less than 5 minutes!**
# 
# *(all timings are for a 1.8GHz Dual-Core Intel Core i5, on a mid-2012 13-inch MacBook Air with 8GB DDR3 memory running macOS Catalina)*
