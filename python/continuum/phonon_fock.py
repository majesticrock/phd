import numpy as np
import matplotlib.pyplot as plt

def plot_function(x_value, y_min=0.9, y_max=1.1, num_points=50000):
    """
    Plots the function:
    sqrt(1+x)*ln|(y+sqrt(1+x))/(y-sqrt(1+x))| - sqrt(1-x)*ln|(y+sqrt(1-x))/(y-sqrt(1-x))|

    Parameters:
    y_value (float): The fixed value for y.
    x_min (float): The minimum value of x.
    x_max (float): The maximum value of x.
    num_points (int): Number of points for plotting.
    """
    # Generate x values within the range
    y = np.linspace(y_min, y_max, num_points)
    
    # Avoid invalid ranges (e.g., values of x that cause division by zero)
    epsilon = 0
    sqrt_y_plus_x = np.sqrt(y + x_value)
    sqrt_y_minus_x = np.sqrt(y - x_value)

    # Compute the function values safely
    try:
        term1 = sqrt_y_plus_x  * np.log(np.abs((1 + sqrt_y_plus_x)  / ((1 - sqrt_y_plus_x)  + epsilon)))
        term2 = sqrt_y_minus_x * np.log(np.abs((1 + sqrt_y_minus_x) / ((1 - sqrt_y_minus_x) + epsilon)))
        f = term1 - term2
    except Exception as e:
        print(f"Error while computing function values: {e}")
        return

    # Plot the function
    plt.figure(figsize=(10, 6))
    plt.plot(y, f, label=f"x = {x_value}")
    plt.title(r"$\sqrt{y+x}\ln\left|\frac{1+\sqrt{y+x}}{1-\sqrt{y+x}}\right| - \sqrt{y-x}\ln\left|\frac{1+\sqrt{y-x}}{1-\sqrt{y-x}}\right|$")
    plt.xlabel("$y$")
    plt.ylabel("$f(y)$")
    #plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    #plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_function_heatmap(x_min=-0.1, x_max=0.1, y_min=0.9, y_max=1.1, num_points=500):
    """
    Plots a heatmap of the function:
    sqrt(1+x)*ln|(y+sqrt(1+x))/(y-sqrt(1+x))| - sqrt(1-x)*ln|(y+sqrt(1-x))/(y-sqrt(1-x))|

    Parameters:
    x_min (float): The minimum value of x.
    x_max (float): The maximum value of x.
    y_min (float): The minimum value of y.
    y_max (float): The maximum value of y.
    num_points (int): Number of points for plotting along each axis.
    """
    # Generate x and y values within the range
    x = np.linspace(x_min, x_max, num_points)
    y = np.linspace(y_min, y_max, num_points)
    X, Y = np.meshgrid(x, y)

    # Avoid invalid ranges (e.g., values of x that cause division by zero)
    epsilon = 1e-8
    sqrt_y_plus_x = np.sqrt(Y + X)
    sqrt_y_minus_x = np.sqrt(Y - X)

    # Compute the function values safely
    try:
        term1 = sqrt_y_plus_x  * np.log(np.abs((1 + sqrt_y_plus_x)  / (1 - sqrt_y_plus_x) + epsilon))
        term2 = sqrt_y_minus_x * np.log(np.abs((1 + sqrt_y_minus_x) / (1 - sqrt_y_minus_x) + epsilon))
        Z = term1 - term2
    except Exception as e:
        print(f"Error while computing function values: {e}")
        return

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    plt.contourf(X, Y, Z, levels=100, cmap='seismic')
    plt.colorbar(label="f(x, y)")
    plt.title(r"$\sqrt{y+x}\ln\left|\frac{1+\sqrt{y+x}}{1-\sqrt{y+x}}\right| - \sqrt{y-x}\ln\left|\frac{1+\sqrt{y-x}}{1-\sqrt{y-x}}\right|$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.grid(False)
    plt.show()


plot_function(x_value=0.002)
#plot_function_heatmap()