'''
Python file to plot 3D energy surfaces as a demonstration of my conformational search workflow for my Honours seminar talk.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from QM.ackley import ackley


def plot_ackley_3d(fig_name):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(-1, 4, 5/2000)
    Y = np.arange(-4, 2, 6/2000)
    X, Y = np.meshgrid(X, Y)

    a = 10
    b = 0.2
    c = 2 * np.pi

    sum_sq_term = -a * np.exp(-b * np.sqrt(X*X + Y*Y) / 2)
    cos_term = -np.exp((np.cos(c*X) + np.cos(c*Y)) / 2)
    Z = a + np.exp(1) + sum_sq_term + cos_term

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.jet, lw=0.5, rstride=50, cstride=50, antialiased=False, alpha=1)
    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.view_init(elev=27, azim=195)
    # fig = plt.gcf()
    # fig.set_size_inches(20, 16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(fig_name, transparent=True)
    plt.show()


def plot_ackley_2d():
    x = np.arange(-32, 33)
    y = [ackley(np.array([x_i])) for x_i in x]
    plt.plot(x,y)
    plt.show()


if __name__ == "__main__":
    plot_ackley_3d("PES_3")
