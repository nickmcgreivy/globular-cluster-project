import csv
from matplotlib import pyplot as plt

energy, time = [] , []

with open('energy_values.csv') as csvfile1:
  read_energy = csv.reader(csvfile1, delimiter = ",")
  
  for row in read_energy:
    energy.append(float(row[7]))
    time.append(float(row[9]))

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.plot(time,energy)

plt.xticks([])
plt.yticks([])
plt.show()


