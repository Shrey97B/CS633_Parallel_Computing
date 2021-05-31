import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

tot_proc = 4
tot_size = 7
tot_exec = 5

procs = [16,36,49,64]
sizes = [16,32,64,128,256,512,1024]

mulsend_data = np.zeros(shape = (tot_proc,tot_size,tot_exec))
pack_data = np.zeros(shape = (tot_proc,tot_size,tot_exec))
derived_data = np.zeros(shape = (tot_proc,tot_size,tot_exec))

#Get timing data from data files and populate them in above numpy arrays
for i in range(tot_proc):
    for j in range(tot_size):
        fil = open("./data_" + str(procs[i]) + "_" + str(sizes[j]) + ".txt","r")
        tim_data = fil.read()
        tim_data = tim_data.split('\n')
        for k in range(tot_exec):
            mulsend_data[i][j][k] = float(tim_data[k*3])
            pack_data[i][j][k] = float(tim_data[k*3+1])
            derived_data[i][j][k] = float(tim_data[k*3+2])

#print(mulsend_data)
#print(pack_data)
#print(derived_data)

#Initialize transparent colors
transp_red = [ 234/256, 162/256, 147/256 ,0.6]
transp_blue = [ 122/256, 237/256, 247/256,0.6 ]
transp_green = [ 155/256, 249/256, 40/256,0.6 ]

for i in range(tot_proc):

    #Form plot for Multiple Send in Blue color
    datap1 = mulsend_data[i]
    datap1 = np.log(datap1)
    med1 = np.median(datap1,axis=1)
    datap1 = np.transpose(datap1)
    #print(datap1)
    plt.title("Plot for Number of Processes = " + str(procs[i]))
    plt.boxplot(datap1,patch_artist=True,boxprops=dict(facecolor=transp_blue, color=transp_blue),
            flierprops=dict(color='blue', markeredgecolor='blue'))
    plt.plot(range(1,tot_size+1),med1,color='blue')
    
    #Form plot for PACK data in Red color
    datap2 = pack_data[i]
    datap2 = np.log(datap2)
    med2 = np.median(datap2,axis=1)
    datap2 = np.transpose(datap2)
    #print(datap2)
    plt.boxplot(datap2,patch_artist=True,boxprops=dict(facecolor=transp_red, color=transp_red),
            flierprops=dict(color='red', markeredgecolor='red'))
    plt.plot(range(1,tot_size+1),med2,color='red')
    
    #Form plot for Derived method data in Green color
    datap3 = derived_data[i]
    datap3 = np.log(datap3)
    med3 = np.median(datap3,axis=1)
    datap3 = np.transpose(datap3)
    #print(datap3)
    plt.boxplot(datap3,patch_artist=True,boxprops=dict(facecolor=transp_green, color=transp_green),
            flierprops=dict(color='green', markeredgecolor='green'))
    plt.plot(range(1,tot_size+1),med3,color='green')
    plt.xlabel('Row/Column Size (N)')
    plt.ylabel('Log(Time)')
    plt.xticks(range(1,tot_size+1),sizes)
    
    #Form legend for plot
    plt.legend(["Multiple Sends", "Using Packed for Data","Using Derived for Data"])
    ax = plt.gca()
    legobj = ax.get_legend()
    legobj.legendHandles[0].set_color('blue')
    legobj.legendHandles[1].set_color('red')
    legobj.legendHandles[2].set_color('green')
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    #Save plot into file with appropriate Name
    fig.savefig('plot' + str(procs[i]) + '.jpg', dpi=100)
    #plt.show()
    plt.close()

print('Plot file plot*.jpg generated')
