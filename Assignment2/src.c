#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>

#define NUM_HOSTS 100

int num_elem,myrank, num_pr, nc_rank, numNode, ppn,host_id;
int nlg_rank,gl_rank;
int nlg_size,gl_size;
int grp_size;
MPI_Comm node_comm, nodel_comm;
MPI_Comm nlg_comm, gl_comm;
long *arr_sizes;

long getRandomValue(long minn,long maxn){
    long x = (long) rand();
    long range = (maxn - minn + 1);
    return (x%range + minn);
}

double getRnDouble(double minn,double maxn){
    double range = maxn - minn;
    double v = ((double)rand()/(double)(RAND_MAX)) * range + minn;
}

void getRandomlyAllocatedSizes(long num_pr, long nume, long* arrs){

    for(int i=0;i<num_pr*num_pr;i++){
        long randomv = getRandomValue(nume/2,nume);
        arrs[i] = randomv;
    }

}

int getHostId(char* hostn,int len){
/*Since there is only one number in Hostid string, it can be extracted through simple concatenation of digits */
    int Hid = 0;
    for(int i=0;i<len;i++){
        if(hostn[i]>='0' && hostn[i]<='9'){
            Hid = Hid*10 + (int)(hostn[i] - '0');
        }
    }
    return Hid;
}

int getGroupId(int hid){
    int grs[6][17] = {{1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,31,-1},
                {13,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,-1},
                {33,34,35,36,37,38,39,40,41,42,43,44,46,-1,-1,-1,-1},
                {45,47,48,49,50,51,52,53,54,56,58,59,60,61,-1,-1,-1},
                {62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78},
                {79,80,81,82,83,84,85,86,87,88,89,90,91,92,-1,-1,-1}};

    for(int i=0;i<6;i++){
        for(int j=0;j<17;j++){
            if(grs[i][j] == hid){
                return (i+1);
            }
        }
    }

    return -1;

}

void MPI_Bcast_default(double* buf,long nume){
    MPI_Bcast(buf,nume,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MPI_Reduce_default(double* inp_buf,long nume,double* recv_buf){
    MPI_Reduce(inp_buf, recv_buf, nume, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void MPI_Gather_default(double* inp_buf,long nume,double* recv_buf){
    MPI_Gather(inp_buf,nume,MPI_DOUBLE,recv_buf,nume,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MPI_Alltoallv_default(double* sendall_def,double* recvall_def,int* sendc,int* recvc,int* sendd,int* recvd){
    MPI_Alltoallv(sendall_def,sendc,sendd,MPI_DOUBLE,recvall_def,recvc,recvd,MPI_DOUBLE,MPI_COMM_WORLD);
}

void MPI_Bcast_optimized(double* Bcast_opt,long nume){

    int ncr = nc_rank;

    if(grp_size==1){
        if(ncr==0)
            MPI_Bcast(Bcast_opt,nume,MPI_DOUBLE,0,nodel_comm);
        if(ppn!=1){
            MPI_Bcast(Bcast_opt,nume,MPI_DOUBLE,0,node_comm);
        }
        return;
    }

    if(nlg_rank==0 && ncr==0){
        MPI_Bcast(Bcast_opt,nume,MPI_DOUBLE,0,gl_comm);
    }
    if(ncr==0 && nlg_size!=1){
        MPI_Bcast(Bcast_opt,nume,MPI_DOUBLE,0,nlg_comm);
    }
    if(ppn!=1){
        MPI_Bcast(Bcast_opt,nume,MPI_DOUBLE,0,node_comm);
    }

}

void MPI_Reduce_optimized(double* Reduce_opt,long nume,double* ReduceRecv_opt){

    double* temp;
    
    if(ppn==1){
        temp = Reduce_opt;
    }
    else{
        temp = (double*) malloc(nume * sizeof(double));
        MPI_Reduce(Reduce_opt,temp, nume, MPI_DOUBLE, MPI_SUM, 0, node_comm);
    }
    if(nc_rank==0){
        MPI_Reduce(temp,ReduceRecv_opt,nume,MPI_DOUBLE,MPI_SUM,0,nodel_comm);
    }

    if(ppn!=1){
        free(temp);
    }

    return;

}

void MPI_Gather_optimized(double* Gather_opt,int nume,double* GatherRecv_opt){
    int ncr=nc_rank;
    double* temp,*temp2;

    if(grp_size==1){
        if(ppn==1){
            temp = Gather_opt;
        }
        else{
            temp = (double*) malloc(nume*ppn*sizeof(double));
            MPI_Gather(Gather_opt,nume,MPI_DOUBLE,temp,nume,MPI_DOUBLE,0,node_comm);
        }

        if(ncr==0){
            MPI_Gather(temp,nume*ppn,MPI_DOUBLE,GatherRecv_opt,nume*ppn,MPI_DOUBLE,0,nodel_comm);
        }

        if(ppn!=1){
            free(temp);
        }

        return;
    }

    //Assuming we are using helper script, every process wil have same num of nodes and every node has same number of process
    int proc_perg = num_pr/grp_size;
    if(ppn==1){
        temp = Gather_opt;
    }
    else{
        temp = (double*) malloc(nume*ppn*sizeof(double));
        MPI_Gather(Gather_opt,nume,MPI_DOUBLE,temp,nume,MPI_DOUBLE,0,node_comm);
    }

    if(ncr==0){
        temp2 = (double*) malloc(nume*proc_perg*sizeof(double));
        MPI_Gather(temp,nume*ppn,MPI_DOUBLE,temp2,nume*ppn,MPI_DOUBLE,0,nlg_comm);
    }

    if(ncr==0 && nlg_rank==0){
        MPI_Gather(temp2,nume*proc_perg,MPI_DOUBLE,GatherRecv_opt,nume*proc_perg,MPI_DOUBLE,0,gl_comm);   
    }

    if(ppn!=1){
        free(temp);
    }
    if(ncr==0){
        free(temp2);
    }


}

void MPI_Alltoallv_optimized(double *sendbuf,double* recvbuf,int* sc,int* rc,int* sd,int* rd){

    if(ppn==1){
        MPI_Alltoallv(sendbuf,sc,sd,MPI_DOUBLE,recvbuf,rc,rd,MPI_DOUBLE,MPI_COMM_WORLD);
        return;
    }

    int *rc_gv,*rd_gv;
    int *sc_gv;
    int *tot_rv;
    int tot_rvt = 0;
    double* gat_out;
    rc_gv = (int*) malloc(num_pr*ppn*sizeof(int));
    rd_gv = (int*) malloc(num_pr*ppn*sizeof(int));
    tot_rv = (int*) malloc(num_pr*sizeof(int));
    if(nc_rank==0){
        for(int i=0;i<num_pr;i++){
            tot_rv[i] = 0;
            for(int j=0;j<ppn;j++){
                rc_gv[i*ppn + j] = arr_sizes[(myrank + j)*num_pr + i];
                rd_gv[i*ppn + j] = tot_rv[i];
                tot_rv[i]+=rc_gv[i*ppn + j];
            }
            tot_rvt+=tot_rv[i];
        }
        gat_out = (double*) malloc(tot_rvt*sizeof(double));
        if(gat_out==NULL){
            printf("Error in memory allocation of gat_out, num elements reqd=%d\n",tot_rvt);
        }
    }

    sc_gv = (int*) malloc(num_pr*sizeof(int));
    int temp_sg = 0;
    int temp_rg = 0;
    for(int i=0;i<num_pr;i++){
        sc_gv[i] = sc[i];
        MPI_Gatherv(sendbuf+temp_sg,sc_gv[i],MPI_DOUBLE,gat_out+temp_rg,rc_gv + (i*ppn),rd_gv+(i*ppn),MPI_DOUBLE,0,node_comm);
        temp_sg+=sc[i];
        temp_rg+=tot_rv[i];
    }

    //if(myrank==0){for(int i=0;i<tot_rvt;i++){printf("%lf ",gat_out[i]);} printf("\n");}

    int *sc_a2a,*sd_a2a;
    int *rc_a2a,*rd_a2a;
    int totsc_a2a = 0;
    int totrc_a2a = 0;
    double* a2a_out;

    if(nc_rank==0){
        sc_a2a = (int*) malloc(numNode*sizeof(int));
        sd_a2a = (int*) malloc(numNode*sizeof(int));
        rc_a2a = (int*) malloc(numNode*sizeof(int));
        rd_a2a = (int*) malloc(numNode*sizeof(int));
        for(int i=0;i<numNode;i++){
            sc_a2a[i] = 0;
            for(int j=0;j<ppn;j++){
                for(int k=0;k<ppn;k++){
                    sc_a2a[i] += rc_gv[i*ppn*ppn + j*ppn + k];
                }
            }
            sd_a2a[i] = totsc_a2a;
            totsc_a2a+=sc_a2a[i];
        }
        for(int i=0;i<numNode;i++){
            rc_a2a[i] = 0;
            for(int j=0;j<ppn;j++){
                for(int k=0;k<ppn;k++){
                    rc_a2a[i] += arr_sizes[(i*ppn + k)*num_pr + (myrank + j)];
                }
            }
            rd_a2a[i] = totrc_a2a;
            totrc_a2a+=rc_a2a[i];
        }
        //if(myrank==0){for(int i=0;i<numNode;i++){printf("%d,%d ",sc_a2a[i],rc_a2a[i]);}printf("\n");}
        a2a_out = (double*) malloc(totrc_a2a*sizeof(double));
        if(a2a_out==NULL){
            printf("Error in memory allocation of a2a_out, num elements reqd=%d\n",totrc_a2a);
        }
        MPI_Alltoallv(gat_out,sc_a2a,sd_a2a,MPI_DOUBLE,a2a_out,rc_a2a,rd_a2a,MPI_DOUBLE,nodel_comm);
        //if(myrank==0){for(int i=0;i<totrc_a2a;i++){printf("%lf ",a2a_out[i]);}printf("\n");}
    }

    int *sc_sca,*sd_sca;
    int rc_sca=0,rd_sca=0;
    int tots_sca = 0;
    int prevs_sca = 0;
    int prevr_sca = 0;

    sc_sca = (int*) malloc(ppn*sizeof(int));
    sd_sca = (int*) malloc(ppn*sizeof(int));

    int ind_sca = 0;
    for(int i=0;i<numNode;i++){
        if(nc_rank==0){
            for(int j=0;j<ppn;j++){
                sc_sca[j] = 0;
                for(int k=0;k<ppn;k++){
                    sc_sca[j]+=arr_sizes[(i*ppn + k)*num_pr + (myrank+j)];
                }

                sd_sca[j] = ((j==0)?0:(sd_sca[j-1]+sc_sca[j-1]));
                tots_sca+=sc_sca[j];
            }
        }
        rc_sca = 0;
        for(int j=0;j<ppn;j++){
            rc_sca += rc[ind_sca];
            ind_sca++;
        }
        MPI_Scatterv(a2a_out+prevs_sca,sc_sca,sd_sca,MPI_DOUBLE,recvbuf+prevr_sca,rc_sca,MPI_DOUBLE,0,node_comm);
        prevr_sca+=rc_sca;
        prevs_sca=tots_sca;
    }

    //if(myrank==0){for(int i=0;i<prevr_sca;i++){printf("%lf ",recvbuf[i]);}printf("\n");}

    if(nc_rank==0){
        free(gat_out);
        free(a2a_out);
    }

}

int main(int argc, char* argv[]){

    //Check command line args

    if(argc!=2){
        printf("Please provide data size as command line argument\n");
        return 0;
    }
    srand(time(NULL));
    int data_size = atoi(argv[1]);
    num_elem = data_size*((int)1024)/sizeof(double);

    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_pr);

//Forming additional communicators for optimization
    double st_opt_ov = MPI_Wtime();

    char* host_n = (char*) malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));
    int hostn_l;
    MPI_Get_processor_name (host_n, &hostn_l);
    host_id = getHostId(host_n,hostn_l);
    MPI_Comm_split(MPI_COMM_WORLD,host_id,myrank,&node_comm);
    MPI_Comm_rank(node_comm,&nc_rank);
    int lead_c = 1;
    if(nc_rank==0){lead_c=0;}
    MPI_Comm_split(MPI_COMM_WORLD,lead_c,myrank,&nodel_comm);
    MPI_Comm_size(node_comm,&ppn);
    numNode = num_pr/ppn;

    int gr_id = getGroupId(host_id);
    MPI_Comm_split(MPI_COMM_WORLD,gr_id*ppn + nc_rank,myrank,&nlg_comm);
    MPI_Comm_size(nlg_comm,&nlg_size);
    lead_c=1;
    MPI_Comm_rank(nlg_comm,&nlg_rank);
    if(nlg_rank==0 && nc_rank==0){
        lead_c=0;
    }

    MPI_Comm_split(MPI_COMM_WORLD,lead_c,myrank,&gl_comm);
    MPI_Comm_rank(gl_comm,&gl_rank);
    MPI_Comm_size(gl_comm,&gl_size);

    if(myrank==0){
        grp_size=gl_size;
    }
    MPI_Bcast(&grp_size,1,MPI_INT,0,MPI_COMM_WORLD);

    double et_opt_ov = MPI_Wtime();
    double dur_opt_ov = et_opt_ov - st_opt_ov, mdur_opt_ov;
    double max_opt;
    MPI_Reduce(&dur_opt_ov,&max_opt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    //printf("Node communication rank for rank %d is rank %d\n",myr,nc_rank);
    if(myrank==0){printf("Time for optimization overhead %lf\n",max_opt);}
    //int numNode = leadc;
    //int ppn = nump/numNode;

    double *Bcast_inp = (double *) malloc(num_elem * sizeof(double));
    for(long i=0;i<num_elem;i++){
        Bcast_inp[i] = getRnDouble(-1000.0,1000.0);
    }

    double *Bcast_def = (double *) malloc(num_elem * sizeof(double));
    if(myrank==0){
        for(long i=0;i<num_elem;i++){
            Bcast_def[i] = Bcast_inp[i];
        }
    }

    double *Bcast_opt = (double *) malloc(num_elem * sizeof(double));
    if(myrank==0){
        for(long i=0;i<num_elem;i++){
            Bcast_opt[i] = Bcast_inp[i];
        }
    }

    double avg1_def = 0.0;
    
    for(int i=0;i<5;i++){
        double st1_def = MPI_Wtime();
        MPI_Bcast_default(Bcast_def,num_elem);
        double et1_def = MPI_Wtime();
        double dur1_def = et1_def - st1_def;
        avg1_def += dur1_def;
    }

/*
    for(long i=0;i<num_elem;i++){
        printf("%lf ",Bcast_def[i] - Bcast_inp[i]);
    }
    printf("\n");*/

    avg1_def = avg1_def/((double) 5.0);
    double mdur1_def;
    MPI_Reduce(&avg1_def,&mdur1_def,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Average Time for Default Broadcast: %lf\n", mdur1_def);
    }

//TODO: code for MPI_Bcast_optimized

    double avg1_opt = 0.0;
    
    for(int i=0;i<5;i++){
        double st1_opt = MPI_Wtime();
        MPI_Bcast_optimized(Bcast_opt,num_elem);
        double et1_opt = MPI_Wtime();
        double dur1_opt = et1_opt - st1_opt;
        avg1_opt += dur1_opt;
    }

/*
    if(myrank==0){
        for(int i=0;i<num_elem;i++){
            printf("%lf ",Bcast_opt[i] - Bcast_def[i]);
        }
        printf("\n");
    }*/

    avg1_opt+=dur_opt_ov;
    avg1_opt = avg1_opt/((double) 5.0);

    double mdur1_opt;
    MPI_Reduce(&avg1_opt,&mdur1_opt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Average Time for Optimized Broadcast: %lf\n", mdur1_opt);
    }

    double *Reduce_inp = (double *) malloc(num_elem * sizeof(double));
    for(long i=0;i<num_elem;i++){
        Reduce_inp[i] = getRnDouble(-1000.0,1000.0);
    }

    double *Reduce_def = (double *) malloc(num_elem * sizeof(double));
    double* ReduceRecv_def;
    for(long i=0;i<num_elem;i++){
        Reduce_def[i] = Reduce_inp[i];
    }

    if(myrank==0){
        ReduceRecv_def = (double *) malloc(num_elem * sizeof(double));
    }

    double *Reduce_opt = (double *) malloc(num_elem * sizeof(double));
    double* ReduceRecv_opt;
    for(long i=0;i<num_elem;i++){
        Reduce_opt[i] = Reduce_inp[i];
    }

    if(myrank==0){
        ReduceRecv_opt = (double *) malloc(num_elem * sizeof(double));
    }

    double avg2_def = 0.0;

    for(int i=0;i<5;i++){
        double st2_def = MPI_Wtime();
        MPI_Reduce_default(Reduce_def,num_elem,ReduceRecv_def);
        double et2_def = MPI_Wtime();
        double dur2_def = et2_def - st2_def;
        avg2_def += dur2_def;
    }

    avg2_def = avg2_def/((double) 5.0);
    double mdur2_def;
    MPI_Reduce(&avg2_def,&mdur2_def,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);


/* if(myrank==0){
    for(long i=0;i<num_elem;i++){
        printf("%lf ",ReduceRecv_def[i]);
    }
    printf("\n");
  }*/


    if(myrank==0){
        printf("Average Time for Default Reduce: %lf\n", mdur2_def);
    }

//TODO: code for MPI_Reduce_optimized
    double avg2_opt = 0.0;

    for(int i=0;i<5;i++){
        double st2_opt = MPI_Wtime();
        MPI_Reduce_optimized(Reduce_opt,num_elem,ReduceRecv_opt);
        double et2_opt = MPI_Wtime();
        double dur2_opt = et2_opt - st2_opt;
        avg2_opt += dur2_opt;
    }

    avg2_opt += dur_opt_ov;
    avg2_opt = avg2_opt/((double) 5.0);
    double mdur2_opt;
    MPI_Reduce(&avg2_opt,&mdur2_opt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

/*
    if(myrank==0){
    for(long i=0;i<num_elem;i++){
        printf("%lf ",ReduceRecv_def[i] - ReduceRecv_opt[i]);
    }
    printf("\n");
    }*/


    if(myrank==0){
        printf("Average Time for Optimized Reduce: %lf\n", mdur2_opt);
    }

    double *Gather_inp = (double *) malloc(num_elem * sizeof(double));
    for(long i=0;i<num_elem;i++){
        Gather_inp[i] = getRnDouble(-1000.0,1000.0);
    }

    double *Gather_def = (double *) malloc(num_elem * sizeof(double));
    double *Gather_opt = (double *) malloc(num_elem * sizeof(double));
    double* GatherRecv_def;
    double* GatherRecv_opt;
    for(long i=0;i<num_elem;i++){
        Gather_def[i] = Gather_inp[i];
    }
    for(long i=0;i<num_elem;i++){
        Gather_opt[i] = Gather_inp[i];
    }

    if(myrank==0){
        GatherRecv_def = (double*) malloc(num_elem * num_pr * sizeof(double));
        GatherRecv_opt = (double*) malloc(num_elem * num_pr * sizeof(double));
    }


    double avg3_def = 0.0;
    for(int i=0;i<5;i++){
        double st3_def = MPI_Wtime();
        MPI_Gather_default(Gather_def,num_elem,GatherRecv_def);
        double et3_def = MPI_Wtime();
        double dur3_def = et3_def - st3_def;
        avg3_def+=dur3_def;
    }


/*   if(myrank==0){
    for(long i=0;i<num_elem*num_pr;i++){
        printf("%lf ",GatherRecv_def[i]);
    }
    printf("\n");
}*/
    avg3_def = avg3_def/((double) 5.0);
    double mdur3_def;
    MPI_Reduce(&avg3_def,&mdur3_def,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Average Time for Default Gather: %lf\n", mdur3_def);
    }

//TODO: code for MPI_Gather_optimized

    double avg3_opt = 0.0;
    for(int i=0;i<5;i++){
        double st3_opt = MPI_Wtime();
        MPI_Gather_optimized(Gather_opt,num_elem,GatherRecv_opt);
        double et3_opt = MPI_Wtime();
        double dur3_opt = et3_opt - st3_opt;
        avg3_opt+=dur3_opt;
    }


/*   if(myrank==0){
    for(long i=0;i<num_elem*num_pr;i++){
        printf("%lf ",GatherRecv_opt[i] - GatherRecv_def[i]);
    }
    printf("\n");
}*/

    avg3_opt+=dur_opt_ov;
    avg3_opt = avg3_opt/((double) 5.0);
    double mdur3_opt;
    MPI_Reduce(&avg3_opt,&mdur3_opt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);    

    if(myrank==0){
        printf("Average Time for Optimized Gather: %lf\n", mdur3_opt);
    }

    arr_sizes = (long *) malloc(num_pr * num_pr* sizeof(long));
    if(myrank==0){

        getRandomlyAllocatedSizes(num_pr,num_elem,arr_sizes);
        /*
        for(int i=0;i<num_pr*num_pr;i++){
            printf("%ld ",arr_sizes[i]);
        }
        printf("\n");*/
    }

    MPI_Bcast(arr_sizes,num_pr*num_pr,MPI_LONG,0,MPI_COMM_WORLD);

    int totscount = 0;
    int totrcount = 0;

    int* send_c = (int*) malloc(num_pr * sizeof(int));
    int* recv_c = (int*) malloc(num_pr * sizeof(int));
    int* send_d = (int*) malloc(num_pr * sizeof(int));
    int* recv_d = (int*) malloc(num_pr * sizeof(int));

    for(int i=0;i<num_pr;i++){
        int s_ind = myrank*num_pr + i;
        send_c[i] = (int)arr_sizes[s_ind];
        send_d[i] = totscount;
        totscount+= (int)arr_sizes[s_ind];

        int r_ind = i*num_pr + myrank;
        recv_c[i] = (int)arr_sizes[r_ind];
        recv_d[i] = totrcount;
        totrcount+= (int)arr_sizes[r_ind];

    }
    /*
    for(long i=0;i<num_pr;i++){
        printf("%ld,%ld,%ld,%ld:%ld ",send_c[i],recv_c[i],send_d[i],recv_d[i],myrank);
    }
    printf("\n");*/

    double* sendall_def = (double*) malloc(totscount * sizeof(double));
    double* recvall_def = (double*) malloc(totrcount * sizeof(double));

    double* sendall_opt = (double*) malloc(totscount * sizeof(double));
    double* recvall_opt = (double*) malloc(totrcount * sizeof(double));

    for(long i=0;i<totscount;i++){
        sendall_def[i] = getRnDouble(-1000.0,1000.0);
    }
    for(int i=0;i<totscount;i++){
        sendall_opt[i] = sendall_def[i];
    }


    double avg4_def = 0.0;

    for(int i=0;i<5;i++){
        if(myrank==0){printf("Default Alltoallv call %d\n",i+1);}
        double st4_def = MPI_Wtime();
        MPI_Alltoallv_default(sendall_def,recvall_def,send_c,recv_c,send_d,recv_d);
        double et4_def = MPI_Wtime();
        double dur4_def = et4_def - st4_def;
        avg4_def+=dur4_def;
    }

    avg4_def = avg4_def/((double) 5.0);
    double mdur4_def;
    MPI_Reduce(&avg4_def,&mdur4_def,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if(myrank==0){
        printf("Average Time for Default Alltoallv: %lf\n", mdur4_def);
        /*
        for(int i=0;i<totscount;i++){
            printf("%lf ",sendall_def[i]);
        }
        printf("\n");
        for(int i=0;i<totrcount;i++){
            printf("%lf ",recvall_def[i]);
        }
        printf("\n");*/
    }

    double avg4_opt = 0.0;

    for(int i=0;i<5;i++){
        if(myrank==0){printf("Optimized Alltoallv call %d\n",i+1);}
        double st4_opt = MPI_Wtime();
        MPI_Alltoallv_optimized(sendall_opt,recvall_opt,send_c,recv_c,send_d,recv_d);
        double et4_opt = MPI_Wtime();
        double dur4_opt = et4_opt - st4_opt;
        avg4_opt+=dur4_opt;
    }

    avg4_opt+=dur_opt_ov;
    avg4_opt = avg4_opt/((double) 5.0);

    double mdur4_opt;
    MPI_Reduce(&avg4_opt,&mdur4_opt,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    //for(int i=0;i<totrcount;i++){printf("%lf ",recvall_opt[i]-recvall_def[i]);}printf("\n");
    if(myrank==0){
        printf("Average Time for Optimized Alltoallv: %lf\n", mdur4_opt);
    }

    if(myrank==0){
        char fileName[100];
        sprintf(fileName,"data_%d_%d.txt",num_pr,data_size);
        FILE *fp = fopen(fileName,"ab+");
        fprintf(fp,"%lf\n%lf\n%lf\n%lf\n",mdur1_def,mdur2_def,mdur3_def,mdur4_def);
        fprintf(fp,"%lf\n%lf\n%lf\n%lf\n",mdur1_opt,mdur2_opt,mdur3_opt,mdur4_opt);
        fclose(fp);

    } 

    MPI_Finalize();
    return 0;
}
