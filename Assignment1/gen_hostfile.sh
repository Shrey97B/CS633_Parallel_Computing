#This is a secondary hostfile generator, which is to be used only when NodeAllocator does not work
#This is a simplistic generator that only checks server liveliness using ping
#To have it work, one will need to uncomment some lines mentioned in run.py

FileName="hosts"
rm -f $FileName
Proc_Per_Node=8
touch $FileName
declare -a NodeArr
Countv=0
echo "Inspecting hosts, please wait..."
for i in $(seq 32)
do
        #Ping the host to know its availability
	ping -c 4 csews$i.cse.iitk.ac.in > /dev/null 2>&1
	res=$?
	if [ $res -eq 0 ]
	then
		#echo "csews$i:$Proc_Per_Node">>$FileName
		#Store Hostname into array if pinged successfully
		NodeArr[$Countv]="csews$i:$Proc_Per_Node"
		Countv=$Countv+1
	fi
	printf "\033[1A"
	printf "\033[K"
	echo "Inspected csews$i"
done
#Randomly shuffle the elements of Node Array
Arr2=( $(echo "${NodeArr[@]}" | sed -r 's/(.[^;]*;)/ \1 /g' | tr " " "\n" | shuf | tr -d " " ) )
#Store the array elements into file
echo ${Arr2[@]} | sed 's/ /\n/g' > $FileName
echo "Hostfile Formed"
