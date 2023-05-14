import sys
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

import plotly.express as px
import plotly.graph_objs as go


def linear_fit(points):
    # construct the matrix A and the vector b
    A = np.column_stack((points[:, 0], np.ones(len(points))))
    b = points[:, 1]
    
    # use the least squares method to find the line parameters
    m, c = np.linalg.lstsq(A, b, rcond=None)[0]
    
    return m, c

#Initial values
# "blobs_expBSA_19_Position_1_SIRT20i_1n.inv.csv"
# "test.dat"

try:
    path_to_data, path_to_save = sys.argv[1].split(' ')
    dir_to_save = path_to_save.split('/')[0]
except:
    raise(Exception("Specify path for data"))    

data_frame = pd.read_csv(path_to_data)

x = data_frame['axis-2'].values
y = data_frame['axis-1'].values
z = data_frame['axis-0'].values

points_3d = np.array([x,y,z]).T
points = np.array([x,z]).T
kx, bx = linear_fit(points)
phi_xz= np.arctan(kx)

fx = lambda x : kx*x + bx

points = np.array([y,z]).T
ky, by = linear_fit(points)
phi_yz= -np.arctan(ky)

fy = lambda x : ky*x + by


U_x = np.array([[1.0, 0.0,             0.0],
                [0.0, np.cos(phi_yz), -np.sin(phi_yz)],
                [0.0, np.sin(phi_yz), np.cos(phi_yz)]])

U_y = np.array([[np.cos(phi_xz), 0.0, np.sin(phi_xz)],
                [0.0,            1.0, 0.0           ],
                [-np.sin(phi_xz),0.0, np.cos(phi_xz)]])

points_3d_new = np.dot(np.dot(U_x,U_y),points_3d.T).T

x_new = points_3d_new.T[0]
y_new = points_3d_new.T[1]
z_new = points_3d_new.T[2]

plt.scatter(y, z)
plt.scatter(y_new, z_new)
plt.plot(y,fy(y), c ='r')

plt.xlabel('y')
plt.ylabel('z')
plt.savefig(dir_to_save+"/xy_py.pdf")

plt.scatter(x, z)
plt.scatter(x_new, z_new)
plt.plot(x,fx(x), c ='r')
plt.xlabel('x')
plt.ylabel('z')

plt.savefig(dir_to_save+"/xz_py.pdf")

# fig = px.scatter_3d(x=x_new,y=y_new,z=z_new)
# fig.update_traces(marker_size = 2)
# fig.show()

# create trace for first scatter plot
trace1 = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=2,
        color='blue',
        opacity=0.5
    ),
    name='Init'
)

# create trace for second scatter plot
trace2 = go.Scatter3d(
    x=x_new,
    y=y_new,
    z=z_new,
    mode='markers',
    marker=dict(
        size=2,
        color='red',
        opacity=0.5
    ),
    name='Rotated'
)

# create layout
layout = go.Layout(
    title='Two Scatter Plots',
    scene=dict(
        xaxis=dict(title='X Axis'),
        yaxis=dict(title='Y Axis'),
        zaxis=dict(title='Z Axis')
    )
)

# create figure and add traces to it
fig = go.Figure(data=[trace1, trace2], layout=layout)

# show the plot
fig.show()

points_to_save = points_3d_new/4.
np.savetxt(path_to_save,points_to_save)    