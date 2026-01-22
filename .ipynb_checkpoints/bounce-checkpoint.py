# -*- coding: utf-8 -*-
"""
Created on ...

@author: ...
Description: ...
    
"""
import matplotlib.pyplot as pyp

# Define functions here
def bounce(height=0.5, v0=4, restitution=1):
    # program
    # Initialize variables
    g = 9.8
    dt = 1
    y = height
    v = v0
    t = 0
    # Create Lists
    y_positions = [y]
    time_values = [t]
    # Simulation loop
    while t < 10:  # Run for 10 seconds
        v = v - g * dt  # Update velocity due to gravity
        y = y + v * dt  # Update position

    # Check for bounce
        if y <= 0:
            y = 0
            v = -v * restitution  # Reverse velocity with energy loss

    # Record data for plotting
        y_positions.append(y)
        time_values.append(t)
        t += dt
    # Plot the results
    pyp.figure(figsize=(8, 6))
    pyp.plot(time_values, y_positions)
    pyp.xlabel('Time (s)')
    pyp.ylabel('Height (m)')
    pyp.title('Bouncing Ball Simulation')
    pyp.grid(True)
    pyp.show()
    return time_values, y_positions


    
#The code below runs whatever is in bounce
if __name__ == '__main__':
    v0 = 0.0   # initial velocity (m/s)
    y0 = 10.0  # initial height (m)
    bounce(y0,v0)