#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define MX_DB (999999.99)

MPI_Comm nodewc,nodelc;
int procpn, numNodes;
int ncw_rnk;

int parseHostID(char* hname,int hlen){
    int h_id = 0;
    for(int i=0;i<hlen;i++){
        if(hname[i]>='0' && hname[i]<='9'){
            h_id = h_id*10 + (int)(hname[i] - '0');
        }
    }
    return h_id;
}

void getNumYears(char* linestr,int* nyears){

    int coln = 0;
    for(int i=0;i<strlen(linestr);i++){
        if(linestr[i]==',' || linestr[i]=='\n'){
            coln = coln+1;
        }
    }
    *nyears = coln - 2;
}

void getColumnYears(char* linestr,int *yearsd,int nyear){

    int coln=0;
    //printf("%d\n",nyear);
    char temp[20];
    strcpy(temp,"");
    for(int i=0;i<strlen(linestr);i++){
        if(linestr[i]==',' || linestr[i]=='\n'){
            //printf("%s\n",temp);
            if(coln>=2){
                yearsd[coln-2] = atoi(temp);
                //printf("%d\n",yearsd[coln-2]);
            }
            coln = coln+1;
            strcpy(temp,"");
        }
        else{
            strncat(temp,&linestr[i],1);
        }
    }
    /*
    for(int i=0;i<nyear;i++){
        printf("Year: %d\n",yearsd[i]);
    }*/
}

int main(int argc,char* argv[]){

    MPI_Init(&argc,&argv);

    int myrank, num_pr;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_pr);

    int *datasca_sc;
    int *datasca_sd;
    float *Latd, *Longtd, *temprdt;
    int numStations = 0;
    int num_years = 0;
    int* year_lab;
    int rceil;

    //Process with rank 0 reads the file
    if(myrank==0){
        //printf("Here\n");
        char* fileName = argv[1];
        FILE* inpFile1 = fopen(fileName,"r");
        if(inpFile1 == NULL){
            printf("Input CSV File not found\n");
        }

        //Considering each cell has max 15 characters and that each row has maximum of 100 years
        char file_lin[2000];

        while(fgets(file_lin,2000,inpFile1)){
            numStations++;
        }
        fclose(inpFile1);
        numStations--;
        printf("Number of stations %d\n",numStations);
        Latd = (float*) malloc(numStations*sizeof(float));
        Longtd = (float*) malloc(numStations*sizeof(float)); 

        FILE* inpFile2 = fopen(fileName,"r");

        int rown = 0;

        while(fgets(file_lin,2000,inpFile2)){
            //printf("%s\n",file_lin);
            int coln = 0;
            if(rown==0){
                getNumYears(file_lin,&num_years);
                year_lab = (int*) malloc(num_years*sizeof(int));

                //Row numbers in array rounded to ceil(numstation/num_pr) * num_pr to use calls without vector variants          
                rceil = (numStations + num_pr - 1)/num_pr;
                rceil = rceil*num_pr;
                temprdt = (float*) malloc(rceil * num_years * sizeof(float));
                getColumnYears(file_lin,year_lab,num_years);
            }
            else{
                char temp[20];
                strcpy(temp,"");
                for(int i=0;i<strlen(file_lin);i++){
                    if(file_lin[i]==',' || file_lin[i]=='\n'){
                        if(coln==0){
                            sscanf(temp, "%f", &Latd[rown-1]);
                        }
                        else if(coln==1){
                            sscanf(temp, "%f", &Longtd[rown-1]);
                        }
                        else{
                            sscanf(temp, "%f", &temprdt[(rown-1)*(num_years) + (coln-2)]);
                        }
                        coln = coln+1;
                        strcpy(temp,"");
                    }
                    else{
                        strncat(temp,&file_lin[i],1);
                    }
                }
            }

            rown = rown + 1;

        }
        /*
        for(int i=(numStations-5);i<numStations;i++){
            for(int j=0;j<num_years;j++){
                printf("%lf ",temprdt[i*num_years + j]);
            }
            printf("\n");
        }
        */
    }

    //Employed Barrier so that the time of other processes waiting for rank 0 to read file is not included
    MPI_Barrier(MPI_COMM_WORLD);

    double sttime = MPI_Wtime();

    char* host_n = (char*) malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));
    int hostn_l;

    //Storing Processor Name using the MPI Call
    MPI_Get_processor_name (host_n, &hostn_l);
    int host_id = parseHostID(host_n,hostn_l);

    //Creating node-wise sub-comm for process in same node
    MPI_Comm_split(MPI_COMM_WORLD,host_id,myrank,&nodewc);
    MPI_Comm_rank(nodewc,&ncw_rnk);
    MPI_Comm_size(nodewc,&procpn);
    numNodes = num_pr/procpn;
    int lead_c = 1;
    if(ncw_rnk==0){lead_c=0;}
    //Creating node-leader based subcommunicator
    MPI_Comm_split(MPI_COMM_WORLD,lead_c,myrank,&nodelc);

    int var_info[2];
    var_info[0] = numStations;
    var_info[1] = num_years;

    //Broadcasting information such as number of stations and Number of Years
    if(ncw_rnk==0){
        MPI_Bcast(var_info,2,MPI_INT,0,nodelc);
    }
    MPI_Bcast(var_info,2,MPI_INT,0,nodewc);

    numStations = var_info[0];
    num_years = var_info[1];

    rceil = (numStations + num_pr - 1)/num_pr;
    rceil = rceil*num_pr;

    //if(myrank==0){printf("%d\n",datasca_rc);for(int i=0;i<datasca_rc;i++){printf("%lf ",recvdata[i]);} printf("\n");}
    //if(myrank==0){printf("%d\n",datasca_rc2);for(int i=0;i<datasca_rc2;i++){printf("%lf ",recvdata[i]);}printf("\n");}

    //Scattering the file data using Node-Based subcommunicators
    int datasca_rc2 = (rceil/numNodes) * num_years,datasca_rc1 =(rceil/num_pr) * num_years;
    float *recvdata2,*recvdata;

    //Scatter the data among Nodes
    if(ncw_rnk==0){
        recvdata2 = (float*) malloc(datasca_rc2 * sizeof(float));
        MPI_Scatter(temprdt,datasca_rc2,MPI_FLOAT,recvdata2,datasca_rc2,MPI_FLOAT,0,nodelc);
    }
    recvdata = (float*) malloc(datasca_rc1*sizeof(float));

    //Scatter the data among process of Nodes
    MPI_Scatter(recvdata2,datasca_rc1,MPI_FLOAT,recvdata,datasca_rc1,MPI_FLOAT,0,nodewc);

    float *ylocal_min = (float*) malloc((num_years+1)*sizeof(float));
    float *yearw_global;
    float loc_min_all;
    float glob_min_all;

    //Initiaizing local minimum values to be maximum possible considering the domain
    for(int j=0;j<num_years;j++){
        ylocal_min[j] = MX_DB;
    }

    //Traversing the domain and computing local minima values for all years
    for(int i=0;i<(rceil/num_pr) && (myrank*(rceil/num_pr) + i)<numStations;i++){
        for(int j=0;j<num_years;j++){
            ylocal_min[j] = (ylocal_min[j]>recvdata[i*num_years + j])?recvdata[i*num_years + j]:ylocal_min[j];
        }
    }

    //Computing the global minima values
    loc_min_all = MX_DB;
    for(int i=0;i<num_years;i++){
        loc_min_all = (loc_min_all>ylocal_min[i])?ylocal_min[i]:loc_min_all;
    }

    //Appending the local overall minima at end
    ylocal_min[num_years] = loc_min_all;
    //if(myrank==0){for(int i=0;i<num_years;i++){printf("%lf ",ylocal_min[i]);}printf("\n");}

    if(myrank==0){
        yearw_global = (float*) malloc((num_years+1)*sizeof(float));
    }

    float *yw_global2;
    if(ncw_rnk==0){
        yw_global2 = (float*) malloc((num_years+1)*sizeof(float));
    }

    //Using reduction to store global year-wise minima and overall minima
    MPI_Reduce(ylocal_min,yw_global2,num_years+1,MPI_FLOAT,MPI_MIN,0,nodewc);
    if(ncw_rnk==0){
        MPI_Reduce(yw_global2,yearw_global,num_years+1,MPI_FLOAT,MPI_MIN,0,nodelc);
    }
    if(myrank==0){
        glob_min_all = yearw_global[num_years];
    }

    //if(myrank==0){for(int i=0;i<num_years;i++){printf("%lf ",yearw_global[i]);}printf("\n");}
    //if(myrank==0){printf("%lf\n",glob_min_all);}

    double ettime = MPI_Wtime();
    double durtime = ettime - sttime;
    double maxtime;
    MPI_Reduce(&durtime,&maxtime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    float maxtf = (float) maxtime;

    //Over-writing the necessary data into file
    if(myrank==0){
        printf("Time taken: %f\n",maxtf);
        char fileName[20] = "output.txt";
        FILE *fp = fopen(fileName,"w+");
        for(int i=0;i<num_years;i++){
            if(i==(num_years-1)){
                fprintf(fp,"%.2f\n",yearw_global[i]);
            }
            else{
                fprintf(fp,"%.2f,",yearw_global[i]);
            }
        }
        fprintf(fp,"%.2f\n",glob_min_all);
        fprintf(fp,"%f\n",maxtf);
    }

    MPI_Finalize();
    return 0;
}
