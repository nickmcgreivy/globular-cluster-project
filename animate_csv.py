import csv;

x = []

bound = 150

with open('test_sim_values.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter = ',')

  for row in readCSV:
    x.append(row)
    print(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

i=0
def animate(i):
  i = 10*i
  if i % 20:
    ax.clear()

  for k in range(len(x[i])//3):
    points = ax.scatter(float(x[i][k]),float(x[i][k+1]),float(x[i][k+2]),color="blue")

  
  ax.set_xlim(-bound,bound)
  ax.set_ylim(-bound,bound)
  ax.set_zlim(-bound,bound)


  return points


ani = animation.FuncAnimation(fig, animate, interval=100, blit=False, frames = 200)
ani.save('orbit.gif', fps=30)

#plt.show()