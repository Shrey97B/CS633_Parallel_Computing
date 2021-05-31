import os

#Remove data files, plot files
os.system("rm -f ./data*.txt")
os.system("rm -f ./plot*.jpg")

#Use Make to remove old object file and compile the src.c file
os.system("make")

numexecs = 10

P_arr = [4,16]
ppn_arr = [1,8]
Ds_arr = [16,256,2048]

G_map = {4:2,16:2}

for exec_num in range(numexecs):
    for P in P_arr:
        for ppn in ppn_arr:
            numG = G_map[P]
            numNG = P//numG
            os.system("python3 script.py " + str(numG) + " " + str(numNG) + " " + str(ppn))
            for D in Ds_arr:
                os.system("mpirun -np " + str(P*ppn) + " -f hostfile ./optim " + str(D))

os.system("python3 create_plots.py")