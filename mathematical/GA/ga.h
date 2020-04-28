#ifndef __GA_H_INCLUDED__
#define __GA_H_INCLUDED__
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include "../test_problem.h"

using namespace std;
# define EPS 1.0e-14

class ga{
	public:
		typedef vector<int> d1;
        typedef vector<d1> d2;
        typedef vector<d2> d3;
        typedef vector<double> dd1;
        typedef vector<dd1> dd2;
		ga(int,int,int,int,double,double,string);
		dd1 run();
	private:
		void init(d3&);
		void crossover(d3&);
        void mutation(d3&);
        void fitness(dd2&);
        void TS(d3&);
        dd2 decode(d3&);
        double evaluation(dd1&);
        

private:
    
	int numRuns;
	int numIter;
	dd1 obj_val;
    dd1 best_sol;
    int population_num;
	double mutation_pro;
    double crossover_pro;
    string func_name;
    int dim; //每一個解的維度
    vector<double> upperbound;
    vector<double> lowerbound;
    d1 nbits;
    test_problem func;
    double best_obj_value;
    
};

ga::ga(int xNumRuns,
int xNumIter,
int xPopulationNum,
int xdim,
double xcrossover_pro,
double xmutation_pro,
string xfunc_name
)
{
    
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
    population_num=xPopulationNum;
    dim=xdim;
    func_name=xfunc_name;
    func=test_problem(func_name,dim);

    upperbound=func.ub;
    lowerbound=func.lb;
    
    crossover_pro=xcrossover_pro;
    
    mutation_pro=xmutation_pro;
    

}
ga::dd2 ga::decode(d3& all_sol){
    dd2 real_sol(population_num,dd1(dim));
    int j, k;
    int sum;
    for(int i=0;i<all_sol.size();i++)
    {
        for (j=0; j<dim; j++)
        {
            sum=0;
            for (k=0; k<nbits[j]; k++)
            {
                if (all_sol[i][j][k]==1)
                {
                    sum += pow(2,nbits[j]-1-k);
                }
            }
           real_sol[i][j] = lowerbound[j] + (double)sum*(upperbound[j] - lowerbound[j])/(double)(pow(2,nbits[j])-1);
        }
    }
    return real_sol;
}

void ga::init(d3& currentSol){
    double rnd;
    nbits.assign(dim,30);
    best_obj_value=numeric_limits<double>::infinity();
    obj_val.assign(population_num,0);
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dim; j++)  
            for (int k=0; k<nbits[j]; k++)
            {
                rnd=(double)rand()/RAND_MAX;
                if (rnd <= 0.5)
                    currentSol[i][j][k] = 0;              
                else               
                    currentSol[i][j][k] = 1;        
            }
}

void ga::TS(d3 &all_sol){
    d3 new_pop(population_num);
    int parent1,parent2;
    //tournament size=2
    for(int i=0;i<population_num;i++){

        parent1=rand()%population_num;
        do{
            parent2=rand()%population_num;
        }while(parent1==parent2);

        if(obj_val[parent1]<obj_val[parent2])
            new_pop[i]=all_sol[parent1];
        else
            new_pop[i]=all_sol[parent2];
    }
    all_sol=new_pop;
}


void ga::crossover(d3 &parent){
    //two-point crossover
    /*
    vector<population>child;
    child=parent;
    double rnd;
    int temp, site1, site2;
    int y=0;
  
    for (int i=0; i<population_num/2; i++)
    {
        for(int j=0;j<dimension;j++){
            rnd = (double)rand()/RAND_MAX;
            if (rnd <= crossover_pro)
            {
                site1 = rand()%(nbits[j]-1);
                site2=rand()%(nbits[j]-1);
                if(site1>site2){
                    temp=site1;
                    site1=site2;
                    site2=temp;
                }
                for(int k=site1;k<=site2;k++){
                    child[y][j][k]=parent[y+1][j][k];
                    child[y+1][j][k]=parent[y][j][k];
                }
    
            }
        }
        y=y+2;
        
    }
    parent=child;
    return;*/
    

    //-------------single point crossover----------------//
    d3 child;
    child=parent;
    double rnd;
    int temp, site1;
    int y=0;
  
    for (int i=0; i<population_num/2; i++)
    {
        for(int j=0;j<dim;j++){
            
            rnd=(double)rand()/RAND_MAX;
            if (rnd <= crossover_pro)
            {
                
                site1=rand()%(nbits[j]-1);
                // site1 = rand()%(nbits[j]-1);
                
                for(int k=0;k<=site1;k++){
                    child[y][j][k]=parent[y+1][j][k];
                    child[y+1][j][k]=parent[y][j][k];
                }
    
            }
        }
        y=y+2;
        
    }
    parent=child;
    return;
    
    //不分維度各維度一起crossover
    // vector<population>child;
    // child=parent;
    // double rnd;
    // int temp, site1;
    // int y=0;
    // int chrom=0,c,w_dim,w_bit;
    // for(int i=0;i<dimension;i++)
    //     chrom+=nbits[i];
    // for (int i=0; i<population_num/2; i++)
    // {

    //     rnd = (double)rand()/RAND_MAX;
    //     if (rnd <= crossover_pro)
    //     {
    //         rnd=(double)rand()/RAND_MAX;
    //         c=floor(rnd*(chrom+10));
            
    //         site1 = c;
    //         if(site1>=chrom)
    //             site1=site1/2;
            
    //         for(int j=0;j< chrom;j++){
                
    //             w_dim=j/nbits[0];           
    //             w_bit=j%nbits[0];
    //             if(j>site1-1){    
    //                 child[y][w_dim][w_bit]=parent[y+1][w_dim][w_bit];        
    //                 child[y+1][w_dim][w_bit]=parent[y][w_dim][w_bit];
    //             }
    //         }
    //     }
    //     y=y+2;  
    // }
    // parent=child;
    // return;
    

    // vector<population>child;
    // child=parent;
    // double rnd;
    // int temp, site1;
    // int y=0;
  
    // for (int i=0; i<population_num/2; i++)
    // {
    //     std::uniform_real_distribution<double> dis(0,1);
    //     rnd=dis(gen);
    //     // rnd = (double)rand()/RAND_MAX;
    //     if (rnd <= crossover_pro)
    //     {
    //         for(int j=0;j<dimension;j++){

    //         // site1 = rand()%(nbits[j]-1);
    //         std::uniform_int_distribution<int> gene(0,nbits[j]-2);
    //         site1=gene(gen);
    //             for(int k=0;k<=site1;k++){
    //                 child[y][j][k]=parent[y+1][j][k];
    //                 child[y+1][j][k]=parent[y][j][k];
    //             }
    //         }
    //     }
        
    //     y=y+2;
        
    // }
    // parent=child;
    // return;
}
void ga::mutation(d3& all_sol){

    double prob;
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dim; j++)
            for (int k=0; k<nbits[j]; k++)
            {
                prob=(double)rand()/RAND_MAX;
                // prob = (double)rand()/RAND_MAX;
                if (prob <=mutation_pro)
                    all_sol[i][j][k]=(all_sol[i][j][k]+1)%2;            
            }

    return;

}
void ga::fitness(dd2& real_sol){
    for(int i=0;i<population_num;i++){
        obj_val[i]=evaluation(real_sol[i]);
        if(obj_val[i]<best_obj_value){
            best_obj_value=obj_val[i];
            best_sol=real_sol[i];
        }
    }
}
ga::dd1 ga::run(){
    
    vector<double>iter_obj_avg(numIter,0.0);
    d3 current_sol(population_num,d2(dim));
    dd2 real_sol(population_num,dd1(dim));
    dd1 each_run_best;
    double each_run_best_obj=numeric_limits<double>::infinity();;
   for(int i=0;i<population_num;i++)
       for(int j=0;j<dim;j++)     
           current_sol[i][j].resize(30);
    
    for(int i=0;i<numRuns;i++){
        
        init(current_sol);
        // cout<<"hello"<<endl;
        real_sol=decode(current_sol);
        
        for(int j=0;j<numIter;j++){
            real_sol=decode(current_sol);
            fitness(real_sol);

            TS(current_sol);
            
            crossover(current_sol);
            
            mutation(current_sol);
    
            
            
            iter_obj_avg[j]+=best_obj_value;
            
        }
        
        if(best_obj_value<each_run_best_obj){
            each_run_best_obj=best_obj_value;
            each_run_best=best_sol;
        }
    }
    
    ofstream output("result/ga/"+func_name+".txt");
    for(int i=0;i<numIter;i++)
        output<<i<<" "<<(double)iter_obj_avg[i]/numRuns<<endl;
    output.close();
    

    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(20)<<(double)iter_obj_avg[i]/numRuns<<endl;
    // cout<<"Evaluation Number:"<<evaluate_num<<endl;
    cout<<"Best objective value:"<<best_obj_value<<endl;
    
    return best_sol;

}



double ga::evaluation(dd1& sol) {
    if(func_name=="Ackley")
        return func.Ackley(sol);
    else if(func_name=="Rastrigin")
        return func.Rastrigin(sol);
    else if(func_name=="Sphere")
        return func.Sphere(sol);
    else if(func_name=="Rosenbrock")
        return func.Rosenbrock(sol);
    else if(func_name=="Michalewicz")
        return func.Michalewicz(sol);
  
}

#endif
       