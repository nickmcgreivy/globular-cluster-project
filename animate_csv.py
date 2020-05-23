import csv;

x = []

bound = 40

with open('test_sim_values.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter = ',')

  for row in readCSV:
    x.append(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

i=0
m = 50
def animate(i):
  ax.clear()

  i = m*i

  for k in range((len(x[i]) - 1)//6):
    points = ax.scatter(float(x[i][6*k]),float(x[i][(6*k)+1]),float(x[i][(6*k)+2]),color="blue")

  
  ax.set_xlim(-bound,bound)
  ax.set_ylim(-bound,bound)
  ax.set_zlim(-bound,bound)

  print("Rendered: " , i)

  return points

print(len(x))
frames = len(x) // m
ani = animation.FuncAnimation(fig, animate, interval=100, blit=False, frames = frames)
ani.save('orbit.gif', fps=60)

