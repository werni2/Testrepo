import math 

def IsEqual(value, reference, eps = 1e-8):
    return abs(value - reference) < 0.5 * eps * abs(value + reference) 

def IsLowerOrEqual(value, reference, eps = 1e-8):
    return value < reference or IsEqual(value, reference, eps)

def IsGreaterOrEqual(value, reference, eps = 1e-8):
    return value > reference or IsEqual(value, reference, eps)

def Threshold(value, border):
    return value if abs(value) > border else 0.0

def bisection(a, b, func, eps=1e-8):

    fa = func(a)
    fb = func(b)
    if fa * fb >= 0:
        raise Exception("You have not assumed right a and b")

    while b - a >= eps:
        # Find middle point
        c = 0.5 * (a + b)
        fc = func(c)
 
        # Check if middle point is root
        if fc == 0.0:
            break
 
        # Decide the side to repeat the steps
        elif fc * fa < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    
    return c

def swap(a, b):
    tmp = a
    a = b
    b = a