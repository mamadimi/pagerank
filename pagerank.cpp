/**
 * @title Parallel and Distributed systems - project 4 - Pagerank
 * @author Mamagiannos Dimitrios (7719)
 * @date March 2015
 */
#include "iostream"
#include <stdlib.h>
#include <fstream>
#include "math.h"
#include "pthread.h"
#include <cstdlib>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>


#define THREADS 8
#define d 0.85


void insert_graph(struct node *nodes,char *grph,unsigned long int N);
void initialize_E_P0(struct posibility *posib,unsigned long int N);
void calculate_posibilities(struct node *nodes,struct posibility *posib,unsigned long int N);
void* get_posib(void *a);
void max_diff(struct posibility *posib,char *pos,unsigned long int N);
void read_pos(double *pos,char *final_pos,unsigned long int N);

using namespace std;

struct timeval startwtime, endwtime;
double seq_time;

typedef struct node{
    unsigned long int index;
    unsigned long int rank_in;      /*For node i: rank_in = sum of nodes that go to i*/
    unsigned long int rank_out;     /*For node i: rank_out = sum of nodes that i goes to.*/
    unsigned long int *in_nodes;    /*For node i: in_nodes[j]=k,contains node k which k->i*/
}node;

typedef struct posibility{
    double posibility;
    double E;
    bool is_steady;                  
}posibility;

typedef struct pthread_args{
    node *nodes;
    posibility *posib;
    unsigned long int N;
    unsigned long int *nodes_row;
    unsigned long int *zero_rank_nodes;
    unsigned long int zero_rank_sum;
    unsigned long int bounds[2];
}pthr_args;

int main(int argc,char **argv){
    
    cout<<"./pagerank N graph.txt posibility.txt\n";
    cout<<"Example: ./pagerank 8846 8846.txt pos8846.txt\n";
    
    if(argc<4) {
        cout<<"Less args\n";
        exit(0);
        
    }
    
    unsigned long int N = atoi(argv[1]);
    char *graph = argv[2];
    char *posibilities = argv[3];
    
    node *nodes = new node[N];
    
    posibility *posib = new posibility[N];
    
    //Insert graph
    insert_graph(nodes,graph,N);
    cout<<"Ok with the graph input\n";
    
    //Initialize E,P0
    initialize_E_P0(posib,N);
    cout<<"Ok with initialize_E_P0\n";
    
    //Calcualte posibilities
    
    calculate_posibilities(nodes,posib,N);
    cout<<"Ok with calculate_posibilities\n";
    
    //Find max_difference of matlab's posibilities and algorithm's posibilities.
    max_diff(posib,posibilities,N);
    
    //Free memory
    for(unsigned long int i=0;i<N;i++){
        delete[] nodes[i].in_nodes;
    }
    
    delete[] nodes;
    delete[] posib;
    
    return 0;
}
/**
 * @param nodes a struct node pointer, information about nodes
 *@param grph graph.txt file 
 *@param N Total nodes
 *@return void
 *
 */
void insert_graph(struct node *nodes,char *grph,unsigned long int N){
    for(unsigned long int i=0;i<N;i++){
        nodes[i].index = i;
        nodes[i].rank_out = 0;
        nodes[i].rank_in = 0;
    }
    
    ifstream p;
    
    //Find the size of in_nodes arrays
    p.open(grph,ios::in);
    if(!p){ 
        cout<<"Cannot open file\n";
        exit(0);
    }
    
    unsigned long int out_node,in_node;
    for(unsigned long int j=0;;){
        
        p>>out_node;
        p>>in_node;
        
        nodes[in_node].rank_in +=1;
        nodes[out_node].rank_out +=1;
        if(p.eof()){	
            break;	
        }    
    }
    p.close();
    
    
    //Allocate memory and fill in_nodes arrays
    
    p.open(grph,ios::in);
    if(!p){ 
        cout<<"Cannot open file\n";   
    }
    
    for(unsigned long int i=0;i<N;i++){
        
        nodes[i].in_nodes = new unsigned long int[nodes[i].rank_in];
        nodes[i].rank_in = 0;
    }
    for(unsigned long int j=0;;){
        p>>out_node;
        
        p>>in_node;
        
        nodes[in_node].in_nodes[nodes[in_node].rank_in] = out_node;
        nodes[in_node].rank_in +=1;
        if(p.eof()){	
            break;
        }    
    }
    p.close();
    
}

/**
 * @param posib a struct posibility pointer, information about posibilies for each node
 * @param N The total nodes
 * @return void
 */

void initialize_E_P0(struct posibility *posib,unsigned long int N){
    
    //Each node has the same posibility and E equals to 1/N;
    for(unsigned long int i=0;i<N;i++){
        posib[i].posibility = (double)(1)/(double)(N);
        posib[i].E = posib[i].posibility;
        posib[i].is_steady = false;
    }
    
}

/**Calculate the new posibility for each node.
 * @param nodes a struct node pointer, information about nodes
 * @param posib a struct posibility pointer, information about posibilies for each node
 * @param n The total nodes
 */

void calculate_posibilities(struct node *nodes,struct posibility *posib,unsigned long int N){
    
    /*Initialize pthreads*/
    
    pthread_t posb_th[THREADS];
    pthread_attr_t attr;
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
    
    
    pthread_args pthr_args[THREADS];
    
    
    /*Find zero rank_out nodes to avoid division by zero. Following matlab's code, for each node with
     *      zero rank_out, all the posibilities added with this node posibility divided by N.
     *      First find out the size of them and then allocate memory and fill the array.
     */
    unsigned long int zero_rank_sum=0;
    for(unsigned long int i=0;i<N;i++){
        if(nodes[i].rank_out==0){
            zero_rank_sum++;
        }
    }
    
    unsigned long int *zero_rank_nodes = new unsigned long int[zero_rank_sum];
    zero_rank_sum = 0;
    
    for(unsigned long int i=0;i<N;i++){
        if(nodes[i].rank_out==0){
            zero_rank_nodes[zero_rank_sum] = i;
            zero_rank_sum++;
        }
    }
    
    for(int i=0;i<THREADS;i++){
        
        pthr_args[i].nodes = nodes;
        pthr_args[i].posib = posib;
        pthr_args[i].N = N;
        pthr_args[i].zero_rank_nodes = zero_rank_nodes;
        pthr_args[i].zero_rank_sum = zero_rank_sum;
    }
    for(int i=0;i<THREADS;i++){
        pthr_args[i].bounds[0] = floor(i*N/THREADS);
        if(i != THREADS-1){
            pthr_args[i].bounds[1] = floor((i+1)*N/THREADS);
        }
        else{
            pthr_args[i].bounds[1] = N;
        }
    }
    
    /**The number of steps to steady state*/
    unsigned long int steps = 0;
    
    /**if true, not in steady state*/
    bool is_balanced;
    
    
    //Start timer
    gettimeofday (&startwtime, NULL);
    do{
        int rc[THREADS];
        for(int i=0;i<THREADS;i++){
            rc[i] = pthread_create(&posb_th[i],&attr,get_posib,(void*)(&pthr_args[i]));
            if(rc[i]!=0){
                cout<<"Creating pthreads error";
                exit(0);
            }
        }
        
        for(int i=0;i<THREADS;i++){
            pthread_join(posb_th[i],NULL);
        }
        //If true continue. if false all nodes are in steady state.
        is_balanced = true;
        
        //Find how much nodes are ready
        unsigned long int ready_nodes = 0;
        for(unsigned long int i=0;i<N;i++){
            if(posib[i].is_steady == true){
                ready_nodes++;
            }
        }
        
        if(ready_nodes>=N){
            is_balanced = false;
        }
        steps++;
        
    }while(is_balanced==true);
    
    //Stop timer
    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
    + endwtime.tv_sec - startwtime.tv_sec);
    cout<<"Total time of posiibilities' calculation = "<<seq_time<<" sec\n";
    cout<<"Steps to steady state ="<<steps<<"\n";
    
    //Find the sum of all possibilities
    double sum_p = 0;
    for(int i=0;i<N;i++){
        sum_p += posib[i].posibility;
    }
    
    
    cout<<"\nSum of all posibilities="<<sum_p<<"\n\n";
    
    delete[] zero_rank_nodes;
}


/**
 * This function is used by p-threads to find out the posibilities
 * @param a void pointer
 * @return void
 */

void* get_posib(void *a){
    
    pthread_args *args = (struct pthread_args*)a;
    
    /*Load struct arguments to local arguments*/
    node *nodes = args->nodes;
    posibility *posib = args->posib;
    unsigned long int N = args->N;
    
    unsigned long int bounds[2];
    bounds[0]=args->bounds[0];
    bounds[1]=args->bounds[1];
    
    double zero_rank_posibility = 0;
    for(unsigned long int j=0;j<args->zero_rank_sum;j++){
        zero_rank_posibility += posib[args->zero_rank_nodes[j]].posibility;
    }
    
    for(unsigned long int i=bounds[0];i<bounds[1];i++){
        
        unsigned long int cur_node = nodes[i].index;
        
        //Add zero rank_out posibilities
        double new_posibility = zero_rank_posibility/(double)N;
        
        for(unsigned long int j=0;j<nodes[cur_node].rank_in;j++){
            new_posibility +=(posib[nodes[cur_node].in_nodes[j]].posibility)/(nodes[nodes[cur_node].in_nodes[j]].rank_out);
        }
        
        new_posibility *= d;
        new_posibility += (1-d)*posib[cur_node].E;
        /*Check if is steady state situation*/
        
        if(fabs(new_posibility-posib[cur_node].posibility)<=1e-4){
            posib[cur_node].is_steady = true;
            posib[cur_node].posibility = new_posibility;
        }
        else{
            posib[cur_node].is_steady = false;
            posib[cur_node].posibility = new_posibility;
            
        }
    }
    
}

/**
 * This function finds the max difference between matlab's posibilities and caclulated posibilities divided by matlab's posibility
 * @param posib a struct posibility pointer, information about posibilies for each node
 * @param pos pos.txt file
 * @param N The total nodes
 * @return void
 */

void max_diff(struct posibility *posib,char *pos,unsigned long int N){
    double *pos_vect = new double[N];
    read_pos(pos_vect,pos,N);
    
    double *quant_diff = new double[N] ;
    double max =0;
    
    for(unsigned long int i=0;i<N;i++){
        quant_diff[i] = ((posib[i].posibility - pos_vect[i])/pos_vect[i]);
        if(fabs(quant_diff[i]>max)) max = quant_diff[i];
        
    }
    cout<<"Max quantity difference = "<<100*max<<"%\n";
    
    delete[] pos_vect;
}

/**
 * This function read posibilities calculated by matlab
 * @param pos pos array
 * @param final_pos pos.txt file
 * @param N The total nodes
 * @return void
 */

void read_pos(double *pos,char *final_pos,unsigned long int N){
    ifstream p;
    p.open(final_pos,ios::in);
    if(!p){
        cout<<"Cannot open file\n";
    }
    double i;
    for(unsigned long j=0;j<N;j++){
        p>>i;
        pos[j] = i;
    }
    p.close();
}
