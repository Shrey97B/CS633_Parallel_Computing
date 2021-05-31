
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt

sns.set()

numexecs = 10
call_strs = ['Bcast','Reduce','Gather','All2allv']
Rec_data = [[[] for i in range(numexecs)] for j in range(4)]


for P in [4, 16]:
    for ppn in [1, 8]:
        for D in [16, 256, 2048]:
            fil = open("./data_" + str(P*ppn) + "_" + str(D) + ".txt","r")
            tim_data = fil.read()
            tim_data = tim_data.split('\n')
            ind=0
            while ind<(numexecs*8):
                d_i = ind//8
                Rec_data[0][d_i].append(float(tim_data[ind]))
                ind=ind+1
                Rec_data[1][d_i].append(float(tim_data[ind]))
                ind=ind+1
                Rec_data[2][d_i].append(float(tim_data[ind]))
                ind=ind+1
                Rec_data[3][d_i].append(float(tim_data[ind]))
                ind=ind+1               

for mn in range(4):
    demo_input_format = pd.DataFrame.from_dict({
        "D": [],
        "P": [],
        "ppn": [],
        "mode": [],  # 1 --> optimized, 0 --> standard
        "time": [],
    })
    ind_dl = [0 for i in range(10)]
    for execution in range(numexecs):
        for P in [4, 16]:
            for ppn in [1, 8]:
                for D in [16, 256, 2048]:
                    dl = Rec_data[mn][execution]
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 1, "time": np.log(dl[ind_dl[execution]+1])
                    }, ignore_index=True)
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 0, "time": np.log(dl[ind_dl[execution]])
                    }, ignore_index=True)
                    ind_dl[execution] = ind_dl[execution] + 2

    demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str, demo_input_format["ppn"])))

    print(demo_input_format)
    #sns.catplot(x="(P, ppn)", y="time", data=demo_input_format, kind="box", col="D", hue="mode")
    g = sns.FacetGrid(data=demo_input_format, col="D")
    g.map_dataframe(sns.boxplot,x="(P, ppn)", y="time", hue="mode")
    g.set_axis_labels("(P,ppn)", "Log(Time)")
    g.add_legend()
    g.fig.set_size_inches(18,8)
    plt.savefig('plot_' + call_strs[mn] + '.jpg')
