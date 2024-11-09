import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#define constants
MASS = 0.145
DIAMETER = 7.4e-2
g = 9.81
RHO = 1.2
C_D = 0.35
A = np.pi*(DIAMETER/2)**2

def rad(theta):
    """
    Converts angles from degrees to radians
    :param theta: an angle in degrees
    :return: the angle theta in radians
    """
    return theta*np.pi/180

def miph_to_mps(x):
    """
    Converts speeds from miles per hour to meters per second
    :param x: the speed in miles per hour
    :return: the speed in meters per second
    """
    return x*1609.34/3600

def ft_to_m(x):
    """
    Converts distances from ft to m
    :param x: the distance in ft
    :return: the distance in m
    """
    return x/3.281

def plot_theoretical(initial_speed, launch_angle, y_0):
    """
    Plots the theoretical curve for a projectile (i.e no air resistance) given an initial speed, launch angle, and y value (x_0 = 0 is assumed)
    :param initial_speed: the initial speed of the projectile in meters per second
    :param launch_angle: the launch angle of the projectile in radians
    :param y_0: the initial y value of the projectile
    :return: nothing (this function is specifically for plotting)
    """
    v_x, v_y = initial_speed * np.cos(launch_angle), initial_speed * np.sin(launch_angle)
    t_final = (v_y+np.sqrt(v_y**2+2*g*y_0))/g
    time = np.linspace(0, t_final, 100)
    x = v_x*time
    y = y_0+v_y*time-g*time**2/2
    plt.plot(x, y)

def plot_fig_2_2(x, y, initial_speed, launch_angle):
    """
    Plots figure 2.2 for a given an array of x and y points and the initial speed and launch angle
    :param x: an array of x coordinates for the solution
    :param y: an array of y coordinates for the solution
    :param initial_speed: the initial speed of the projectile in meters per second
    :param launch_angle: the launch angle of the projectile in radians
    :return:
    """
    # code for plotting fig 2.2 (also set y to 0 in solve)
    plot_theoretical(initial_speed, launch_angle, x[-1])
    plt.xlabel("Range (m)")
    plt.ylabel("Height (m)")
    plt.title("Projectile motion")
    plt.legend(["Theoretical (no air resistance)", "Euler Method"])
    plt.xlim(0, 25)
    plt.ylim(-1, 7)
    plt.show()

def plot_fig_2_3(x, y, initial_speed, launch_angle, legend):
    """
    plot figure 2.3 for a given an array of x and y points and the initial speed and launch angle
    :param x: an array of x coordinates for the solution
    :param y: an array of y coordinates for the solution
    :param initial_speed: the initial speed of the projectile in meters per second
    :param launch_angle: the launch angle of the projectile in radians
    :param legend: the text of the legend entry for the ode solving method
    :return:
    """
    plot_theoretical(initial_speed, launch_angle, y[0])
    plt.scatter(x, y, marker="P")
    plt.xlabel("Range (m)")
    plt.ylabel("Height (m)")
    plt.title("Projectile motion")
    plt.legend(["theoretical (no air resistance)", legend])
    plt.xlim(0, 250)
    plt.ylim(-10, 70)
    plt.show()

def solve(initial_speed, launch_angle, time_step, method, air_resistance=True):
    """
    Solves for the motion of a projectile with or without air resistance, and returns the range of the projectile
    :param initial_speed: the initial speed of the projectile in meters per second
    :param launch_angle: the launch angle of the projectile in degrees
    :param time_step: the difference in time between each point in the solution in seconds
    :param method: the solution method to use (euler, euler-cromer, or midpoint)
    :param air_resistance: whether to include air resistance in the solution (True or False)
    :return: the height of the ball (in meters) 400ft from where it was launched, or -10000 if the projectile was ever at a height lower than -10000m
    """
    fence_distance = ft_to_m(400)
    method = method.lower()
    launch_angle = rad(launch_angle)
    lower_y_limit = -10000

    #extract the components of the initial velocity and set the initial position
    x_array, y_array = [], []
    x, y = 0, 1
    v_x, v_y = initial_speed * np.cos(launch_angle), initial_speed * np.sin(launch_angle)

    #condition should be modified to just y>=0 for pt 2
    while x <= fence_distance and y > lower_y_limit:
        x_array.append(x)
        y_array.append(y)

        if air_resistance:
            magv = np.sqrt(v_x**2 + v_y**2)
            a = ((-v_x*C_D*RHO*A*magv)/(2*MASS),
                (-v_y*C_D*RHO*A*magv)/(2*MASS) - g)
        else:
            a = (0, -g)

        if method == "euler":
            x = x + time_step*v_x
            y = y + time_step*v_y
            v_x = v_x + time_step*a[0]
            v_y = v_y + time_step*a[1]
        elif method == "euler-cromer":
            v_x = v_x + time_step*a[0]
            v_y = v_y + time_step*a[1]
            x = x + time_step*v_x
            y = y + time_step*v_y
        elif method == "midpoint":
            v_nx = v_x
            v_ny = v_y
            v_x = v_x + time_step*a[0]
            v_y = v_y + time_step*a[1]
            x = x + time_step*((v_x+v_nx)/2)
            y = y + time_step*((v_y+v_ny)/2)
        else: break

    x_array.append(x)
    y_array.append(y)

    #uncomment to plot figures (Part 1)
    #plot_fig_2_2(x_array, y_array, initial_speed, launch_angle)
    #plot_fig_2_3(x_array, y_array, initial_speed, launch_angle, f"{method} method")

    #calculate range (Part 2)
    #return -y_array[-2]*((x_array[-1]-x_array[-2])/(y_array[-1]-y_array[-2])) + x_array[-2]

    #calculate range (Part 3)
    if y <= lower_y_limit:
        return lower_y_limit
    return y_array[-2]+(((y_array[-1]-y_array[-2])/(x_array[-1]-x_array[-2]))*(fence_distance-x_array[-2]))



#function call for plotting fig 2.2
#solve(15, 45, 0.1, "euler", air_resistance=False)

#function call for plotting fig 2.3 (change method parameter to euler, euler-cromer or midpoint)
#solve(50, 45, 0.1, "midpoint", air_resistance=True)

###Part 2
#generate distributions
samp_num = 10000
time_step = 0.1
velocity_dist = miph_to_mps(15 * np.random.randn(samp_num) + 100)
theta_dist = 10*np.random.randn(samp_num)+45

#compute the AB/HR as specified in part 2
"""ranges = np.empty(samp_num)
for i in range(samp_num):
    ranges[i] = solve(velocity_dist[i], theta_dist[i], time_step, method="euler", air_resistance=True)

AB_HR = samp_num/len(np.where(ranges>ft_to_m(400))[0])
print(f"AB/HR ratio for a fence at 0m: {AB_HR}")"""

###Part 3
fence_height_num = 50
fence_heights = np.linspace(0.5, 15, fence_height_num)
ABHRs = np.empty(fence_height_num)

ranges = np.empty(samp_num)
for i in range(samp_num):
   ranges[i] = solve(velocity_dist[i], theta_dist[i], time_step, method="euler", air_resistance=True)

for i in range(fence_height_num):
    ABHRs[i] = samp_num / len(np.where(ranges > fence_heights[i])[0])

abhr_func = sp.interpolate.interp1d(fence_heights, ABHRs-10)
fence_height = sp.optimize.bisect(abhr_func, 0.5, 15)
print(f"Fence height required for AB/HR <= 10: {fence_height}")
plt.scatter(fence_heights, ABHRs)
plt.title("ABHR vs Fence Height for 10000 baseballs")
plt.xlabel("Fence Height (m)")
plt.ylabel("ABHR")
plt.show()

