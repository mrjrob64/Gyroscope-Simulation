import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

data = np.loadtxt('simple_gyro_v2.csv', delimiter=',')
time_data = data[:, 0]
x_data = data[:, 1]
y_data = data[:, 2]
z_data = data[:, 3]
angular_velocity_x = data[:, 4]
angular_velocity_y = data[:, 5]
angular_velocity_z = data[:, 6]
torque_x = data[:, 7]
torque_y = data[:, 8]
torque_z = data[:, 9]   
velocity_x = data[:, 10]
velocity_y = data[:, 11]
velocity_z = data[:, 12]

STEPS = 100
START = 0

# Function to update the plot for each frame
def update(frame):
    # Clear previous plot
    ax.cla()

    # Plot a single point in 3D space
    ax.scatter([point[0]], [point[1]], [point[2]], c='red', s=50, label='Moving Point')

    # Plot angular velocity vector
    ax.quiver([0], [0], [0], [angular_velocity_x[frame * STEPS + START]], [angular_velocity_y[frame * STEPS + START]], [angular_velocity_z[frame * STEPS + START]], color='green', label='Angular Velocity Vector', length=1, normalize=True)

    # Plot torque vector
    ax.quiver([0], [0], [0], [torque_x[frame * STEPS + START]], [torque_y[frame * STEPS + START]], [torque_z[frame * STEPS + START]], normalize=True, color='orange', label='Torque Vector')

    # plot position vector and gravity vector
    ax.quiver([point[0]], [point[1]], [point[2]], [0], [0], [-1], color='red', label='Gravity Vector', length=0.4, normalize=True)
    ax.quiver([0], [0], [0], [point[0]], [point[1]], [point[2]], color='black', label='Position Vector', length=0.4, normalize=True)

    # Plot line from origin to point
    ax.plot([0, point[0]], [0, point[1]], [0, point[2]], c='black', label='Line from Origin to Point')
    
    #plot velocity vector
    ax.quiver([point[0]], [point[1]], [point[2]], [velocity_x[frame * STEPS + START]], [velocity_y[frame * STEPS + START]], [velocity_z[frame * STEPS + START]], color='orange', label='Velocity', length=0.4)

    # Plot the path
    path_x.append(point[0])
    path_y.append(point[1])
    path_z.append(point[2])
    ax.plot(path_x, path_y, path_z, c='blue', label='Path')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set fixed axis limits
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(0, 1.1)

    # Update the position of the point
    point[0] = x_data[frame * STEPS + START]
    point[1] = y_data[frame * STEPS + START]
    point[2] = z_data[frame * STEPS + START]

# Create a figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial position of the point
point = [x_data[START], y_data[START], z_data[START]]

# Lists to store the path
path_x = [point[0]]
path_y = [point[1]]
path_z = [point[2]]

# Create an animation
animation = FuncAnimation(fig, update, frames=200, interval=1)

# Display the animation
animation.save('simple_gyro.gif', writer='imagemagick', fps=10)
