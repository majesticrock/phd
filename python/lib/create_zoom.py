# ax                    -   axes to plot on, needs to be an object; no plt wrapper
# inset_x/ypos, 
#   inset_width/height  -   position and size of the inset in relative coordinates (0 to 1)
# x/ylim                -   x and y limit of the zoomed region, expects tuples (lower limit, upper limit)
# xy_data               -   tuple of the (x, y)-data. If it is None the data is automatically extracted from the existing lines
# x/yticklabels         -   Default behaviour: Do not show labels in the zoom; can be changed however by specifying them
# **kwargs              -   is forwarded to the creation of the inset axes
def create_zoom(ax, inset_xpos, inset_ypos, inset_width, inset_height, xlim=(0, .1), ylim=(0, .1), xy_data=None, xticklabels=[], yticklabels=[], **kwargs):
    axins = ax.inset_axes([inset_xpos, inset_ypos, inset_width, inset_height], xlim=xlim, ylim=ylim, xticklabels=xticklabels, yticklabels=yticklabels, **kwargs)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    
    for line in ax.lines:
        if xy_data is None:
            # Extract the data from the line
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            new_line, = axins.plot(x_data, y_data)
        else:
            new_line, = axins.plot(xy_data[0], xy_data[1])
        
        # Create a new object for the inset plot
        # I do not know how this works, but it does
        new_line.update_from(line)
        new_line.remove()
        new_line.axes = axins
        new_line.set_transform(axins.transData)
        axins.add_line(new_line)
        
## Usage example
#import matplotlib.pyplot as plt
#import numpy as np
#x = np.linspace(0, 10, 100)
#fig, ax = plt.subplots()
#
#ax.plot(x, np.sin(x),           linestyle="--", label="sin", color="orange", marker="v", markevery=4)
#ax.plot(x, np.cos(x),           linestyle="--", label="cos", color="purple", marker="o", markevery=3)
#ax.plot(x, np.sin(np.sin(x)),   linestyle="--", label="sin(sin)", color="green", marker="X", markevery=5)
#create_zoom(ax, inset_xpos=0.45, inset_ypos=0.25, inset_width=0.25, inset_height=0.35, xlim=(0.65, 0.95), ylim=(0.6, 0.8))
#
#ax.legend()
#ax.set_xlabel("$x$")
#ax.set_ylabel("$y$")
#
#fig.tight_layout()
#plt.show()