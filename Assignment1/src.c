#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

//Function to generate a random double value within min and max
double range_random(double min,double max){
    double ans = (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
    return ans;
}

//Returns the top neighbour of process with rank proc_rank and nump processes
//-1 if such neighbour does not exist
int top_n(int proc_rank,int nump){

    int rowlen = sqrt(nump);
    if(proc_rank - rowlen<0){
        return -1;
    }
    return proc_rank - rowlen;
}

//Returns the bottom neighbour of process with rank proc_rank and nump processes
//-1 if such neighbour does not exist
int bottom_n(int proc_rank,int nump){

    int rowlen = sqrt(nump);
    if(proc_rank + rowlen>=nump){
        return -1;
    }
    return proc_rank + rowlen;
}

//Returns the left neighbour of process with rank proc_rank and nump processes
//-1 if such neighbour does not exist
int left_n(int proc_rank,int nump){

    int rowlen = sqrt(nump);
    if(proc_rank%rowlen==0){
        return -1;
    }
    return proc_rank - 1;
}

//Returns the right neighbour of process with rank proc_rank and nump processes
//-1 if such neighbour does not exist
int right_n(int proc_rank,int nump){

    int rowlen = sqrt(nump);
    if((proc_rank+1)%rowlen==0){
        return -1;
    }
    return proc_rank + 1;
}

//function for copying 2D array
void copy2DArray(double** src,double** targ,int m,int n){
    for(int i=0;i<m;i++){
        memcpy(targ[i],src[i],n*sizeof(double));
    }
}

//check if the given indice values are valid for nxn matrix
int check_valid_index(int i,int j,int n){
    if(i>=0 && i<n && j>=0 && j<n){
        return 1;
    }
    return 0;
}

//Computes Stencil for Pack Method
void compute_Stencil_Pack(double** arr,int pr_rank,int nump,int subd_n,double* buf[4]){

    //create temporary 2D array with current matrix data
    double** temp;
    temp = (double**) malloc(sizeof(double*)*subd_n);
    for(int i=0;i<subd_n;i++){
        temp[i] = (double*) malloc(sizeof(double)*subd_n);
    }

    copy2DArray(arr,temp,subd_n,subd_n);
    
    //pos variable for all 4 neighbours
    int pos[4] = {0,0,0,0};
    
    //rank of 4 neighbours
    int neighb[4];
    neighb[0] = top_n(pr_rank,nump);
    neighb[1] = bottom_n(pr_rank,nump);
    neighb[2] = left_n(pr_rank,nump);
    neighb[3] = right_n(pr_rank,nump);
    
    //below indi,indj values helps provide indices for neighbouring cells
    //i.e. temp[i-1][j],temp[i+1][j],temp[i][j-1],temp[i][j+1]
    int indi[] = {-1,1,0,0};
    int indj[] = {0,0,-1,1};

    for(int i=0;i<subd_n;i++){
        for(int j=0;j<subd_n;j++){
            double currs = 0;
            double currc = 4;

            for(int x=0;x<4;x++){
                if(check_valid_index(i+indi[x],j+indj[x],subd_n)){
                    //The neighbouring cells are present in the process matrix itself
                    currs+=temp[i+indi[x]][j+indj[x]];
                }
                else if(neighb[x]!=-1){
                    //Neighbour cells extracted from neighbouring buffers using Unpack
                    double tn =0;
                    MPI_Unpack(buf[x],sizeof(double),&pos[x],&tn,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    //printf("%lf %d %d\n",tn,i,j);
                    /*
                        if(pr_rank==4){
                            //printf("%lf %d %d\n",tn,i,j);
                        }
                    */
                    currs+=tn;
                }
                else{
                    //Neighbouring cell does not exist due to cell being present at overall domain boundary
                    currc--;
                }
            }
            arr[i][j] = currs/currc;

        }
    }

}


//function computing Stencil for both Multiple Send and Derived method
void compute_Stencil_Mus_Derived(double** arr,int pr_rank,int nump,int subd_n,double* buf[4]){

    //create temporary 2D array with current matrix data
    double** temp;
    temp = (double**) malloc(sizeof(double*)*subd_n);
    for(int i=0;i<subd_n;i++){
        temp[i] = (double*) malloc(sizeof(double)*subd_n);
        for(int j=0;j<subd_n;j++){
            temp[i][j] = arr[i][j];
        }
    }

    //pos variable for all 4 neighbours    
    int pos[4] = {0,0,0,0};
    
    //rank of 4 neighbours
    int neighb[4];
    neighb[0] = top_n(pr_rank,nump);
    neighb[1] = bottom_n(pr_rank,nump);
    neighb[2] = left_n(pr_rank,nump);
    neighb[3] = right_n(pr_rank,nump);
    
    //below indi,indj values helps provide indices for neighbouring cells
    //i.e. temp[i-1][j],temp[i+1][j],temp[i][j-1],temp[i][j+1]
    int indi[] = {-1,1,0,0};
    int indj[] = {0,0,-1,1};

    for(int i=0;i<subd_n;i++){
        for(int j=0;j<subd_n;j++){
            double currs = 0;
            double currc = 4;

            for(int x=0;x<4;x++){
                if(check_valid_index(i+indi[x],j+indj[x],subd_n)){
                    //The neighbouring cells are present in the process matrix itself
                    currs+=temp[i+indi[x]][j+indj[x]];
                }
                else if(neighb[x]!=-1){
                    //Get values from neighbouring buffers received
                    double tn =0;
                    tn = buf[x][pos[x]];
                    pos[x]++;
                    //printf("%lf %d %d",tn,i,j);
                    currs+=tn;
                }
                else{
                    //Cell is on boundary of entire domain
                    currc--;
                }
            }
            arr[i][j] = currs/currc;

        }
    }

}

//Function for communication for Multiple Send
void compute_communicate_mus(double** arr,int t,int myrank,int nump,int subd_n){

    //Send and Receive requests
    MPI_Request *s_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4*subd_n);
    MPI_Request *r_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4*subd_n);
    
    //Status array
    MPI_Status *stat = (MPI_Status*) malloc(sizeof(MPI_Status)*8*subd_n);

    //Send and Receive buffers
    double *buf_s[4];
    double *buf_r[4];
    
    //count to store number of send/receive requests
    int reqs = 0;

    //rank of neighbouring process
    int neighb_ne[4];
    neighb_ne[0] = top_n(myrank,nump);
    neighb_ne[1] = bottom_n(myrank,nump);
    neighb_ne[2] = left_n(myrank,nump);
    neighb_ne[3] = right_n(myrank,nump);

    // below variables for indexing boundary elements to send
    // indexing [indr + posr*i][indc+posc*i] -> [0][i],[subd_n-1][i],[i][0],[i][subd_n-1]
    int indr[] = {0,subd_n-1,0,0};
    int indc[] = {0,0,0,subd_n-1};

    int posr[] = {0,0,1,1};
    int posc[] = {1,1,0,0};

    //iterating over four boundaries
    for(int x=0;x<4;x++){

        if(neighb_ne[x]!=-1){
            buf_s[x] = (double*) malloc(sizeof(double)*subd_n);
            buf_r[x] = (double*) malloc(sizeof(double)*subd_n);
            
            //iterating over each boundary to index and send values
            for(int i=0;i<subd_n;i++){
                buf_s[x][i] = arr[indr[x] + posr[x]*i][indc[x] + posc[x]*i];
                MPI_Isend(&buf_s[x][i],1,MPI_DOUBLE,neighb_ne[x],t*subd_n + i,MPI_COMM_WORLD,&s_req[reqs]);
                MPI_Irecv(&buf_r[x][i],1,MPI_DOUBLE,neighb_ne[x],t*subd_n + i,MPI_COMM_WORLD,&r_req[reqs]);
                reqs++;
            }
        }

    }

    MPI_Waitall(reqs,r_req,stat);
    MPI_Waitall(reqs,s_req,stat+reqs);

    compute_Stencil_Mus_Derived(arr,myrank,nump,subd_n,buf_r);
}

//function for processing PACK communication
void compute_communicate_pack(double** arr,int t,int myrank,int nump,int subd_n){

    //Send and Receive requests
    MPI_Request *s_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4);
    MPI_Request *r_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4);
    MPI_Status *stat = (MPI_Status*) malloc(sizeof(MPI_Status)*8);

    //Buffers for send/receive
    double *buf_s[4];
    double *buf_r[4];
    int reqs = 0;

    int top_ne = top_n(myrank,nump);
    if(top_ne!=-1){
        reqs++;
        buf_s[0] = (double*) malloc(sizeof(double)*subd_n);
        buf_r[0] = (double*) malloc(sizeof(double)*subd_n);
        int pos = 0;
        //Packing the entire top row
        MPI_Pack(arr[0],subd_n,MPI_DOUBLE,buf_s[0],subd_n*sizeof(double),&pos,MPI_COMM_WORLD);
        MPI_Isend(buf_s[0],pos,MPI_PACKED,top_ne,t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[0],subd_n,MPI_DOUBLE,top_ne,t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int bottom_ne = bottom_n(myrank,nump);
    if(bottom_ne!=-1){
        reqs++;
        buf_s[1] = (double*) malloc(sizeof(double)*subd_n);
        buf_r[1] = (double*) malloc(sizeof(double)*subd_n);
        int pos = 0;
        //Packing the bottom row
        MPI_Pack(arr[subd_n-1],subd_n,MPI_DOUBLE,buf_s[1],subd_n*sizeof(double),&pos,MPI_COMM_WORLD);
        MPI_Isend(buf_s[1],pos,MPI_PACKED,bottom_ne,t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[1],subd_n,MPI_DOUBLE,bottom_ne,t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int left_ne = left_n(myrank,nump);
    if(left_ne!=-1){
        reqs++;
        buf_s[2] = (double*) malloc(sizeof(double)*subd_n);
        buf_r[2] = (double*) malloc(sizeof(double)*subd_n);
        int pos = 0;
        //Iterating over left boundary and PACKING the values into buf_s[2]
        for(int i=0;i<subd_n;i++){
            MPI_Pack(&arr[i][0],1,MPI_DOUBLE,buf_s[2],subd_n*sizeof(double),&pos,MPI_COMM_WORLD);
        }
        MPI_Isend(buf_s[2],pos,MPI_PACKED,left_ne,t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[2],subd_n,MPI_DOUBLE,left_ne,t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int right_ne = right_n(myrank,nump);
    if(right_ne!=-1){
        reqs++;
        buf_s[3] = (double*) malloc(sizeof(double)*subd_n);
        buf_r[3] = (double*) malloc(sizeof(double)*subd_n);
        int pos = 0;
        //Iterating over right boundary and PACKING the values into buf_s[3]
        for(int i=0;i<subd_n;i++){
            MPI_Pack(&arr[i][subd_n-1],1,MPI_DOUBLE,buf_s[3],subd_n*sizeof(double),&pos,MPI_COMM_WORLD);
        }
        MPI_Isend(buf_s[3],pos,MPI_PACKED,right_ne,t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[3],subd_n,MPI_DOUBLE,right_ne,t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }


    //Need to wait for receive requests for further computation
    MPI_Waitall(reqs,r_req,stat);


    compute_Stencil_Pack(arr,myrank,nump,subd_n,buf_r);

    //As send buffers are stored separately, can wait for those requests later
    MPI_Waitall(reqs,s_req,stat+reqs);

}

void compute_communicate_Derived(double** arr,int tag,int t,int myrank,int nump,int subd_n,MPI_Datatype rowm,MPI_Datatype colm){

    //Send and Receive requests
    MPI_Request *s_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4);
    MPI_Request *r_req = (MPI_Request*) malloc(sizeof(MPI_Request)*4);
    MPI_Status *stat = (MPI_Status*) malloc(sizeof(MPI_Status)*8);

    //Receive buffers
    double *buf_r[4];
    int reqs = 0;

    int top_ne = top_n(myrank,nump);
    if(top_ne!=-1){
        reqs++;
        buf_r[0] = (double*) malloc(sizeof(double)*subd_n);
        //Sending a rowm type element i.e. (entire row) starting with arr[0]
        MPI_Isend(arr[0],1,rowm,top_ne,tag+t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[0],subd_n,MPI_DOUBLE,top_ne,tag+t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int bottom_ne = bottom_n(myrank,nump);
    if(bottom_ne!=-1){
        reqs++;
        buf_r[1] = (double*) malloc(sizeof(double)*subd_n);
        //Sending a rowm type element i.e. (entire row) starting with arr[N-1]
        MPI_Isend(arr[subd_n-1],1,rowm,bottom_ne,tag+t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[1],subd_n,MPI_DOUBLE,bottom_ne,tag+t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int left_ne = left_n(myrank,nump);
    if(left_ne!=-1){
        reqs++;
        buf_r[2] = (double*) malloc(sizeof(double)*subd_n);
        
	//Sending a colm type element i.e. (entire column) starting with arr[0][0]
        MPI_Isend(&arr[0][0],1,colm,left_ne,tag+t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[2],subd_n,MPI_DOUBLE,left_ne,tag+t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    int right_ne = right_n(myrank,nump);
    if(right_ne!=-1){
        reqs++;
        buf_r[3] = (double*) malloc(sizeof(double)*subd_n);
	//Sending a colm type element i.e. (entire column) starting with arr[0][N-1]
        MPI_Isend(&arr[0][subd_n-1],1,colm,right_ne,tag+t,MPI_COMM_WORLD,&s_req[reqs-1]);
        MPI_Irecv(buf_r[3],subd_n,MPI_DOUBLE,right_ne,tag+t,MPI_COMM_WORLD,&r_req[reqs-1]);
    }

    MPI_Waitall(reqs,r_req,stat);
    //Need to wait for send requests as well since using the same matrix address onto which computation is applied
    MPI_Waitall(reqs,s_req,stat+reqs);
    compute_Stencil_Mus_Derived(arr,myrank,nump,subd_n,buf_r);

}

int main(int argc,char* argv[]){

    if(argc!=3){
        printf("Please provide data points per process and number of time step as arguments\n");
        return 0;
    }

    //Data points per process through cmd arg
    int dp_per_process = atoi(argv[1]);
    //printf("%d",dp_per_process);

    //Number of timestamps (50) through cmd args
    int num_ts = atoi(argv[2]);
    //printf("%d",num_ts);

    int myrank = 0;
    int nump = 0;

    MPI_Init(&argc,&argv);
    
    //process rank in myrank
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
    //number of processes in nump
    MPI_Comm_size(MPI_COMM_WORLD,&nump);

    //Using Square root to get row/col length for matrices
    int subd_n = sqrt(dp_per_process);
    
    //Initializing and populating the input data for process
    double** subdomain;
    subdomain = (double**) malloc(sizeof(double*)*subd_n);
    for(int i=0;i<subd_n;i++){
        subdomain[i] = (double*) malloc(sizeof(double)*subd_n);
    }

    for(int i=0;i<subd_n;i++){
        for(int j=0;j<subd_n;j++){
            subdomain[i][j] = range_random(-5000.0,5000.0);
        }
    }

/*
    printf("%d\n",top_n(myrank,nump));
    printf("%d\n",bottom_n(myrank,nump));
    printf("%d\n",left_n(myrank,nump));
    printf("%d\n",right_n(myrank,nump));*/

    //Using Multiple Send method
    //Variables for recording time
    double stt1,ent1,dur1,maxt1;
    double** array1 = (double**) malloc(sizeof(double*)*subd_n);
    for(int i=0;i<subd_n;i++){
        array1[i] = (double*) malloc(sizeof(double)*subd_n);
    }
    //Copying input values
    copy2DArray(subdomain,array1,subd_n,subd_n);
    stt1 = MPI_Wtime();
    //Iterating Number of Timestamp times
    for(int t=0;t<num_ts;t++){
    //call Multiple send communication method
        compute_communicate_mus(array1,t,myrank,nump,subd_n);
        /*
        for(int i=0;i<subd_n;i++){
            for(int j=0;j<subd_n;j++){
                printf("%lf %d %d %d\n",array1[i][j],myrank,i,j);
            }
        }
        */
    }

    ent1 = MPI_Wtime();
    dur1 = ent1 - stt1;

    //Using reduce to get the maximum time spent by all processes
    MPI_Reduce(&dur1,&maxt1,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Time taken by Multiple Sends: %lf\n",maxt1);
    }

    //Using Pack Method
    
    //Variables for recording time
    double stt2,ent2,dur2,maxt2;

    double** array2;
    array2 = (double**) malloc(sizeof(double*)*subd_n);
    for(int i=0;i<subd_n;i++){
        array2[i] = (double*) malloc(sizeof(double)*subd_n);
    }

    //Copy Input
    copy2DArray(subdomain,array2,subd_n,subd_n);

    stt2 = MPI_Wtime();
    int tag2 = num_ts*subd_n;
    for(int t = 0; t<num_ts;t++){
        //Call Pack communication method

        compute_communicate_pack(array2,tag2 + t,myrank,nump,subd_n);
        

        /*
        for(int i=0;i<subd_n;i++){
            for(int j=0;j<subd_n;j++){
                printf("%lf %d %d %d\n",array2[i][j],myrank,i,j);
            }
        }
        */
    }


    ent2 = MPI_Wtime();
    dur2 = ent2 - stt2;

    //Call reduce to get maximum time duration among all processes
    MPI_Reduce(&dur2,&maxt2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Time taken by MPI Pack: %lf\n",maxt2);
    }


    //Using contiguous to form a matrix row datatype
    MPI_Datatype rowm;
    MPI_Type_contiguous(subd_n,MPI_DOUBLE,&rowm);

    //Using vector to form a matrix column datatype
    MPI_Datatype colm;
    MPI_Type_vector(subd_n,1,subd_n,MPI_DOUBLE,&colm);

    MPI_Type_commit(&rowm);
    MPI_Type_commit(&colm);

    //Variables for recording time
    double stt3,ent3,dur3,maxt3;

    double **array3 = (double**)malloc(subd_n*sizeof(double*));
    array3[0] = (double*)malloc(subd_n*subd_n*sizeof(double));
    for(int ind_i=1;ind_i<subd_n;ind_i++){
        array3[ind_i] = array3[ind_i-1] + subd_n;
    }
    
    //Copying input into array
    copy2DArray(subdomain,array3,subd_n,subd_n);


    //Stencil Computation using Derived Types
    stt3 = MPI_Wtime();
    //Tag values to continue after the above method communication
    int tag3 = num_ts*(subd_n + 1);
    for(int t = 0; t<num_ts;t++){

        //call derived method communication
        compute_communicate_Derived(array3,tag3,t,myrank,nump,subd_n,rowm,colm);

        /*
        for(int i=0;i<subd_n;i++){
            for(int j=0;j<subd_n;j++){
                printf("%lf %d %d %d\n",array3[i][j],myrank,i,j);
            }
        }
        */
    }


    ent3 = MPI_Wtime();
    dur3 = ent3 - stt3;

    //Using Reduce to compute maximum duration among all processes
    MPI_Reduce(&dur3,&maxt3,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Time taken by MPI Derived: %lf\n",maxt3);
    }    

    /*
    for(int i=0;i<subd_n;i++){
        for(int j=0;j<subd_n;j++){
            printf("%lf ",array1[i][j] - array3[i][j]);
        }
        printf("\n");
    }*/

    /*
    for(int i=0;i<subd_n;i++){
        for(int j=0;j<subd_n;j++){
            printf("%lf ",array1[i][j] - array2[i][j]);
        }
        printf("\n");
    }*/

    //Write data into data_nump_N.txt files
    //The data will be appended to old data
    if(myrank==0){
        char fileName[100];
        sprintf(fileName,"data_%d_%d.txt",nump,subd_n);

        FILE *fp = fopen(fileName,"ab+");
        fprintf(fp,"%lf\n%lf\n%lf\n",maxt1,maxt2,maxt3);
        fclose(fp);

    }

    MPI_Finalize();



    return 0;
}
