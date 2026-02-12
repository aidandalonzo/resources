# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 16:42:19 2018

@author: Tom K
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def euler(func, a=0., b=10., N=100):
    """Computes solution using euler method
    
    Parameters
    ----------
    func  : lambda function
        first derivative expression
    a,b  : 
        interval of solution
    N : float
        number of steps
    """
    h = (b-a)/N
    # initial condition
    x = 0
    tpoints = np.arange(a,b,h)
    xpoints = []
    for t in tpoints:
        xpoints.append(x)
        x += h*func(x,t)
    #plots solution
    plt.plot(tpoints,xpoints)
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.show()
    
    return tpoints, xpoints
    

def rk2(func,a=0.,b=10.,N=100):
    """Computes solution using 2nd order runge kutta method
    
    Parameters
    ----------
    func  : lambda function
        first derivative expression
    a,b  : 
        interval of solution
    N : float
        number of steps
    """
    h = (b-a)/N

    tpoints = np.arange(a,b,h)
    xpoints = []
    #sets initial condition
    x = 0.0
    for t in tpoints:
        xpoints.append(x)
        k1 = h*func(x,t)
        k2 = h*func(x+0.5*k1,t+0.5*h)
        x += k2
    #plots solution
    plt.plot(tpoints,xpoints)
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.show()
    return tpoints, xpoints

def rk4(func,a=0.,b=10.,N=100):
    """Computes solution using 4th order runge kutta method
    
    Parameters
    ----------
    func  : lambda function
        first derivative expression
    a,b  : 
        interval of solution
    N : float
        number of steps
    """
    h = (b-a)/N

    tpoints = np.arange(a,b,h)
    xpoints = []
    #sets initial condition
    x = 0.0
    for t in tpoints:
        xpoints.append(x)
        k1 = h*func(x,t)
        k2 = h*func(x+0.5*k1,t+0.5*h)
        k3 = h*func(x+0.5*k2,t+0.5*h)
        k4 = h*func(x+k3,t+h)
        x += (k1+2*k2+2*k3+k4)/6
    #plots solution
    plt.plot(tpoints,xpoints)
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.show()
    return tpoints, xpoints

def radau2(func, a=0., b=10., N=100):
    """1-stage Radau IIA method (implicit, 3rd order, L-stable)
    
    Solves dx/dt = f(x,t) using implicit evaluation at right endpoint.
    This is the simplest Radau IIA method.
    """
    h = (b - a) / N
    tpoints = np.linspace(a, b, N)
    xpoints = []
    
    x = 0.0  
    xpoints.append(x)

    for t in tpoints[:-1]:
        # Implicit equation: x_new = x + h*f(x_new, t+h)
        def implicit_eq(x_new):
            return x_new - x - h * func(x_new, t + h)

        # Solve for x_new
        x_new = fsolve(implicit_eq, x)[0]  # Initial guess: x
        xpoints.append(x_new)
        x = x_new

    return np.array(tpoints), np.array(xpoints)


def bdf2(func, a=0., b=10., N=100):
    """Computes solution using BDF-2 (Backward Differentiation Formula)
    
    Parameters
    ----------
    func  : lambda function
        First derivative expression
    a, b  : float
        Interval of solution
    N : int
        Number of steps
    """
    h = (b - a) / N

    # Time points
    tpoints = np.linspace(a, b, N)
    xpoints = []
    
    # Initial condition
    x = 0.0  
    xpoints.append(x)

    # Use Euler for first step
    x1 = x + h * func(x, tpoints[0])
    xpoints.append(x1)

    # Implicit BDF-2 method for subsequent steps
    for i in range(1, len(tpoints) - 1):
        t = tpoints[i + 1]
        
        # Define the implicit equation to solve
        def implicit_eq(x_new):
            return x_new - (4/3) * xpoints[-1] + (1/3) * xpoints[-2] - (2/3) * h * func(x_new, t)

        # Solve for x_new using fsolve
        x_new = fsolve(implicit_eq, xpoints[-1])[0]
        xpoints.append(x_new)

    # Plot solution
    plt.plot(tpoints, xpoints, label="BDF-2")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.legend()
    plt.grid()
    plt.show()

    return tpoints, xpoints

# The LSODA method is much more involved, this is my best guess on it
def lsoda_approx(func, a=0., b=10., N=100, tol=1e-3):
    """Approximates LSODA by switching between Adams-Bashforth and BDF

    Parameters
    ----------
    func  : lambda function
        First derivative expression
    a, b  : float
        Interval of solution
    N : int
        Number of steps
    tol : float
        Tolerance for stiffness detection
    """
    h = (b - a) / N
    tpoints = np.linspace(a, b, N)
    xpoints = []

    # Initial condition
    x = 0.0
    xpoints.append(x)

    # Use explicit Euler (Adams-Bashforth 1-step) for first step
    x1 = x + h * func(x, tpoints[0])
    xpoints.append(x1)

    # Loop through all time points
    for i in range(1, len(tpoints) - 1):
        t = tpoints[i + 1]
        x_old = xpoints[-1]

        # Estimate stiffness using derivative magnitude
        stiffness_factor = abs(func(x_old, t) - func(xpoints[-2], tpoints[i])) / h

        if stiffness_factor < tol:
            # Use explicit Adams-Bashforth (non-stiff)
            x_new = x_old + h * (3/2 * func(x_old, t) - 1/2 * func(xpoints[-2], tpoints[i]))
        else:
            # Use implicit BDF-2 (stiff)
            def implicit_eq(x_new):
                return x_new - (4/3) * x_old + (1/3) * xpoints[-2] - (2/3) * h * func(x_new, t)

            x_new = fsolve(implicit_eq, x_old)[0]

        xpoints.append(x_new)

    # Plot solution
    plt.plot(tpoints, xpoints, label="LSODA Approximation")
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.legend()
    plt.grid()
    plt.show()

    return tpoints, xpoints



# Define ODE: dx/dt = -2x + cos(t)
f = lambda x, t: -2*x + np.cos(t)

# Compare methods
t1, x1 = euler(f, a=0., b=5., N=100)
t2, x2 = rk2(f, a=0., b=5., N=100)
t3, x3 = rk4(f, a=0., b=5., N=100)
t4, x4 = radau2(f, a=0., b=5., N=100)
t5, x5 = bdf2(f, a=0., b=5., N=100)
t6, x6 = lsoda_approx(f, a=0., b=5., N=100, tol=1e-3)