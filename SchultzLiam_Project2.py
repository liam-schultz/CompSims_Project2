import numpy as np
import matplotlib.pyplot as plt

#define constants
MASS = 0.145
DIAMETER = 7.4e-2
g = 9.81
RHO = 1.2
C_D = 0.35
A = np.pi*(DIAMETER/2)**2

def rad(theta):
    return theta*np.pi/180

def integrate(inital_speed, launch_angle, time_step, method, air_resistance=True):
    method = method.lower()
    launch_angle = rad(launch_angle)

    #extract the components of the initial velocity and set the initial position
    x_array, y_array = [], []
    x, y = 0, 0
    v_x, v_y = inital_speed*np.cos(launch_angle), inital_speed*np.sin(launch_angle)
    while y >= 0:
        x_array.append(x)
        y_array.append(y)
        if air_resistance:
            a = ((-abs(v_x)*C_D*RHO*A*abs(v_x))/(2*MASS),\
                (-abs(v_y)*C_D*RHO*A*abs(v_y))/(2*MASS) - g)
        else:
            a = (0, -g)

        if method == "euler":
            x = x+time_step*v_x
            y = y+time_step*v_y
            v_x = v_x+time_step*a[0]
            v_y = v_y+time_step*a[1]
    x_array.append(x)
    y_array.append(y)

    #code for plotting fig 2.2 (also set y to 0 above)
    plt.scatter(x_array, y_array, marker="P")
    plt.xlabel("Range (m)")
    plt.ylabel("Height (m)")
    plt.xlim(0, 25)
    plt.ylim(-1, 7)
    plt.show()

integrate(15, 45, 0.1, "euler", air_resistance=False)

