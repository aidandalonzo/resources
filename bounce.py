# -*- coding: utf-8 -*-
"""
Created on 1/19/26
@author: aidandalonzo
Description: Bouncing Ball 2D
"""

import matplotlib.pyplot as pyp

def bounce(height, v0, restitution):
    # Constants
    g = 9.8
    dt = 0.1

    # Initialize variables
    x = 0
    y = height
    vx = v0[0]
    vy = v0[1]
    t = 0
    
    stop_height = 0.1 * height
    last_peak = height
    current_peak = 0
    going_up = False

    # Create Lists
    x_positions = [x]
    y_positions = [y]
    time_values = [t]

    # Simulation loop
    while t < 100:  # Run for 100 seconds
        vy = vy - g * dt  # Update velocity due to gravity

        # Update positions
        x = x + vx * dt
        y = y + vy * dt
        t += dt

        # Detect upward motion
        if vy > 0:
            going_up = True
        if y > current_peak:
            current_peak = y
        
        # Check for bounce
        if y <= 0:
            y = 0
            vy = -vy * restitution # Reverse velocity with energy loss
            last_peak = current_peak
            current_peak = 0
            going_up = False

        # Stop when peak is too small
        if last_peak < stop_height:
            break

        # Record data for plotting
        x_positions.append(x)
        y_positions.append(y)
        time_values.append(t)

    # Plot y vs time
    pyp.figure(figsize=(8, 6))
    pyp.plot(time_values, y_positions)
    pyp.xlabel("Time (s)")
    pyp.ylabel("Height (m)")
    pyp.title("Y Position vs Time")
    pyp.grid(True)

    # Plot y vs x
    pyp.figure(figsize=(8, 6))
    pyp.plot(x_positions, y_positions)
    pyp.xlabel("X Position (m)")
    pyp.ylabel("Y Position (m)")
    pyp.title("Y Position vs X Position")
    pyp.grid(True)

    pyp.show()

    return time_values, x_positions, y_positions


# Run the program
if __name__ == "__main__":
    v0 = [4.0, 0.0]   # [vx, vy]
    y0 = 10.0
    restitution = 0.8
    bounce(y0, v0, restitution)
