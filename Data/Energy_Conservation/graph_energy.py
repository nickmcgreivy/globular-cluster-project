import csv
from matplotlib import pyplot as plt

energy, time = [],[]

with open('Data/Energy_Conservation/Energy.csv') as csvfile:
  read_energy = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy:
    energy.append(float(row[7]))
    time.append(float(row[9]))

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.plot(time,energy)

ax1.set_title("Leap Frog with Adaptive Time-Stepping Algorithm Energy Conservation")
ax1.set_xlabel("Time")
ax1.set_ylabel("Change in Total Energy of System")
plt.show()


