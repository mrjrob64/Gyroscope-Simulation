import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt('simple_gyro.csv', delimiter=',')
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

STEPS = 1
START = 500000

# Function to update the plot for each frame
def update(frame, point, path_x, path_y):
    # Clear previous plot
    plt.clf()

    # Plot a single point in 2D space
    plt.scatter(x_data[frame * STEPS + START], y_data[frame * STEPS + START], c='red', s=50, label='Moving Point')

    # Plot angular velocity vector in 2D
    plt.quiver(x_data[frame * STEPS + START], y_data[frame * STEPS + START], angular_velocity_x[frame * STEPS + START], angular_velocity_y[frame * STEPS], color='green', label='Angular Velocity Vector')

    # Plot torque vector in 2D
    plt.quiver(x_data[frame * STEPS + START], y_data[frame * STEPS + START], torque_x[frame * STEPS + START], torque_y[frame * STEPS + START], color='orange', label='Torque Vector')

    # plot position vector and gravity vector in 2D
    plt.quiver(x_data[frame * STEPS + START], y_data[frame * STEPS + START], 0, 0, color='red', label='Gravity Vector', scale=0.4)
    plt.quiver(0, 0, x_data[frame * STEPS + START], y_data[frame * STEPS + START], color='black', label='Position Vector', scale=0.4)

    # Plot the path in 2D
    path_x.append(x_data[frame * STEPS + START])
    path_y.append(y_data[frame * STEPS + START])
    plt.plot(path_x, path_y, c='blue', label='Path')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim(-0.11, 0.11)
    plt.ylim(-0.11, 0.11)

    # Display legend
    plt.legend()

# Create a figure
fig = plt.figure()

# Initial position of the point
point = [x_data[0], y_data[0]]

# Lists to store the path
path_x = [point[0]]
path_y = [point[1]]

# Create an animation
animation = FuncAnimation(fig, update, frames=500, fargs=(point, path_x, path_y), interval=1)

# Display the animation
animation.save('simple_gyro_2d.gif', writer='imagemagick', fps=10)

plt.show()
