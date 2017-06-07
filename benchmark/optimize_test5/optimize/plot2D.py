import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from matplotlib import cm, ticker
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

df = pd.read_csv("gaoptimization (copy).csv")
x = df["Phase_Tolerance"].as_matrix()
y = df["Angular_Tolerance"].as_matrix()
z = df["distance"].as_matrix()

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

fig = plt.figure(1)
ax = fig.add_subplot(121, projection='3d')
surf = ax.plot_surface(xi,yi,zi, cmap = plt.cm.jet,antialiased=True)
ax.set_xlabel("Phase Tolerance")
ax.set_ylabel("Angular Tolerance")
ax.set_zlabel("Distance")

ax = fig.add_subplot(122)
#pcolormesh
plt.contourf(xi, yi, zi, cmap = plt.cm.jet)
ax.set_xlabel("Phase Tolerance")
ax.set_ylabel("Angular Tolerance")
plt.colorbar()
plt.show()
