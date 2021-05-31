import os

#Remove data files, plot files
os.system("rm -f ./data*.txt")
os.system("rm -f ./plot*.jpg")

#Use Make to remove old object file and compile the src.c file
os.system("make")

Parr = [16,36,49,64]
Narr = [16**2,32**2,64**2,128**2,256**2,512**2,1024**2]

#compile the Node Allocator
os.system("g++ $HOME/UGP/allocator/src/allocator_improved.cpp -o $HOME/UGP/allocator/src/allocator.out")

for execution in range(1,6):
    #Use the Node Allocator .obj file to create host file
    os.system("$HOME/UGP/allocator/src/allocator.out 64 8")
    for P in Parr:
        for N in Narr:
            #Uncomment below two lines to use a simpler hostfile generator, only if NodeAllocator does not work
            #os.system("chmod 700 ./gen_hostfile.sh")
            #os.system("./gen_hostfile.sh")
            os.system("mpirun -np " + str(P) + " -f hosts ./halo " + str(N) + " 50")

#Execute the plots creator script
os.system("python3 create_dplots.py")
