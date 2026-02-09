# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:07:02 2025

@author: tomke and claude
"""

import numpy as np
import matplotlib.pyplot as plt

def forward(y, x):
    """
    Computes the forward difference derivative.
    """
    dy = np.diff(y)  
    dx = np.diff(x)  
    return dy / dx, x[:-1]  

def backward(y, x):
    """
    Computes the backward difference derivative.
    )
    """
    dy = y[1:] - y[:-1] 
    dx = x[1:] - x[:-1] 
    return dy / dx, x[1:]  

def central(y, x):
    """
    Computes the central difference derivative.
    """
    dy = y[2:] - y[:-2]  
    dx = x[2:] - x[:-2]  
    return dy / dx, x[1:-1]  

if __name__ == '__main__':
    # Test script
    x = np.linspace(0, 10, 100)
    y = np.sin(x)  # Function to differentiate
    
    fwd, x_fwd = forward(y, x)
    bwd, x_bwd = backward(y, x)
    ctr, x_ctr = central(y, x)
    
    plt.figure(figsize=(8, 5))
    plt.plot(x_fwd, fwd, label='Forward Difference', linestyle='dashed')
    plt.plot(x_bwd, bwd, label='Backward Difference', linestyle='dotted')
    plt.plot(x_ctr, ctr, label='Central Difference', linestyle='solid')
    plt.plot(x, np.cos(x), label='Analytical Derivative', linestyle='dashdot')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel("Derivative")
    plt.title("Numerical Differentiation")
    plt.grid()
    plt.show()
