import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Removing current output files
os.system("rm -f ./output*.txt")

#Creating new overall output file for storing output of every config
os.system("touch output_all.txt")

#Compiling the code
os.system("make")

numNodes = [1,2]
corePN = [1,2,4]

#Fixing the input csv file name
DataFileName = 'tdata.csv'

plot_tick = []

#Edit below variable to change number of executions
numexecs = 5

for execs in range(numexecs):
    for nn in numNodes:
        for cpn in corePN:
            print('Execution number:', str(execs+1) + ', Configuration: Number of nodes=', str(nn), ', Process per Node=', str(cpn))
            print('Generating Hostfile')
            os.system("python3 hostf_gen.py 1 " + str(nn) + " " + str(cpn))
            os.system("mpirun -np " + str(nn*cpn) + " -f hostfile ./code " + DataFileName)
            os.system("cat output.txt >> output_all.txt")

for nn in numNodes:
    for cpn in corePN:   
        plot_tick.append("(" + str(nn) + "," + str(cpn) + ")") 

#Fetching output time data
time_list = [[0 for j in range(numexecs)] for i in range(len(numNodes)*len(corePN))]
fil = open("./output_all.txt","r")
file_data = fil.read()
file_data = file_data.split('\n')
#print(file_data)
for exs in range(numexecs):
    for i in range(len(numNodes)*len(corePN)):
        time_list[i][exs] = float(file_data[3*(exs*len(numNodes)*len(corePN) + i) + 2])
print(time_list)

#Plotting line chart for time-readings
tnp = np.array(time_list)
tnp_med = np.median(tnp,axis=1)
plt.boxplot(time_list)
plt.plot(range(1,len(plot_tick)+1),tnp_med)
plt.xticks(range(1,len(plot_tick)+1),plot_tick)
plt.xlabel('Num Nodes, PPN configuration')
plt.ylabel('Time')
plt.savefig('plot.png')

