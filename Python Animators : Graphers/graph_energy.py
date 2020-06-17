import csv
from matplotlib import pyplot as plt

energy_LF, time_LF, energy_AT, time_AT = [],[], [], []

with open('Data/Output/EnergyLF.csv') as csvfile:
  read_energy_LF = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy_LF:
    energy_LF.append(float(row[7]))
    time_LF.append(float(row[9]))

with open('Data/Output/EnergyAT.csv') as csvfile:
  read_energy_AT = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy_AT:
    energy_AT.append(float(row[7]))
    time_AT.append(float(row[9]))

fig = plt.figure()
ax1 = fig.add_subplot(111)

line1, = ax1.plot(time_LF,energy_LF,color="blue")
line2, = ax1.plot(time_AT,energy_AT,color="red")
	
line1.set_label("RK")
line2.set_label("Adaptive Stepping")

ax1.legend()

ax1.set_title("Leap Frog with Adaptive Time-Stepping Algorithm Energy Conservation")
ax1.set_xlabel("Time")
ax1.set_ylabel("Change in Total Energy of System")
plt.show()


