from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



def first_derivative(x,y):
	delx = x[1]-x[0]
	dy_list = np.copy(y)

	for i in range(len(y)):

		if i==0:
			dy = 0

		elif i == len(y) - 1:
			dy = 0

		else:

			dy = (y[i+1] - y[i-1])/(2*delx)

		dy_list[i] = dy

	return dy_list

def second_derivative(x,y):
	delx = x[1]-x[0]
	dy_list = np.copy(y)

	for i in range(len(y)):

		if i==0:
			dy = (y[1] - y[0]) / delx

		elif i == len(y) - 1:
			dy = (y[i] - y[i-1]) / delx

		else:

			dy = (y[i+1] - y[i-1])/(2*delx)

		dy_list[i] = dy

	return dy_list


def local_laplacian(temps, x, y):
	dtemps = np.copy(temps)

	dx = np.copy(x)
	for i in range(len(x)):
		dx[i] = second_derivative(x[i],first_derivative(x[i],temps[i]))
		
	dy = np.copy(y)	
	for i in range(len(y)):
		dy[:,i] = second_derivative(y[:,i],first_derivative(y[:,i],temps[:,i]))


	for i in range(len(y)):
		for k in range(len(x)):
			dtemps[i][k] = dx[i][k] + dy[i][k]

	return dtemps
		

def animate(i):
	global z, two
	
	dt = 0.13

	ax.clear()
	
	z = z + heat_diffusion*local_laplacian(z, x, y)*dt
	color = cm.jet(z)
	surface = ax.plot_surface(x,y,z, facecolors = color)
			

	ax.set_xlim([x_limits[0],x_limits[1]])
	ax.set_ylim([y_limits[0],y_limits[1]])
	ax.set_zlim([z_limits[0],z_limits[1]])

	print("Successfully rendered frame: ", i)

	return surface


x_limits = [-2,2]
y_limits = [-2,2]
z_limits = [-3,3]
resolution = 40
heat_diffusion = 0.005

fig = plt.figure()
ax = fig.add_subplot(111, projection ='3d')


ax.set_xlim([x_limits[0],x_limits[1]])
ax.set_ylim([y_limits[0],y_limits[1]])
ax.set_zlim([z_limits[0],z_limits[1]])

ax.set_zlim([-2,4])
x = np.linspace(x_limits[0],x_limits[1], resolution)
y = np.linspace(y_limits[0],y_limits[1], resolution)

x,y = np.meshgrid(x,y)

two = np.copy(x)
for i in range(len(x)):
	for k in range(len(x)):
		two[i][k] = 1

#z = 4/(two+(x**2 + y**2)**2)
z = np.sin(x**2 + y**2)

ani = animation.FuncAnimation(fig, animate, interval=100, blit=False, frames = 400)

#ani.save('heat_diffusion2.gif', fps=45)

plt.show()






