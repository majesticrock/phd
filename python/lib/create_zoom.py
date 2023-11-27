def create_zoom(ax, inset_xpos, inset_ypos, inset_width, inset_height, xlim=(0, .1), ylim=(0, .1), xticklabels=[], yticklabels=[], **kwargs):
    axins = ax.inset_axes([inset_xpos, inset_ypos, inset_width, inset_height], xlim=xlim, ylim=ylim, xticklabels=xticklabels, yticklabels=yticklabels, **kwargs)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    
    for line in ax.lines:
        # Extract the data from the line
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        
        # Create a new object for the inset plot
        # I do not know how this works, but it does
        new_line, = axins.plot(x_data, y_data)
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