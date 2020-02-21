#ifndef __GOA_H_INCLUDED__
#define __GOA_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include<iomanip>
using namespace std;
#define pi 3.141592653;

class goa{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		goa(int,int,string,int,double,double,int);
		solution run();
	private:
		void init();
        double dis_cal(solution,solution);
        void generate_new_solutions();
		double evaluation(solution);
        double frand();
        double deg_forces(double);
        solution transition(population,int,double);
        double coef_c(int);

        void fitness(population&);

private:
	int numRuns;
	int numIter;
    int dim;
	string filename;
	population currentSol;
	solution all_obj_value;
	double best_obj_value;
	solution best_sol;
    int population_num;
    double upperbound;
    double lowerbound;
    double cmax;
    double cmin;
};

goa::goa(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
double xcmax,
double xcmin,
int xdim
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	filename=xfilename;
    population_num=xPopulationNum;
    cmax=xcmax;
    cmin=xcmin;
    dim=xdim;
    upperbound=32.768;
    lowerbound=-32.768;
    
}

double goa::deg_forces(double r){
    double f=0.5;
    double l=1.5;
    return f*exp((-r)/l)-exp(-r)+0.0000184043;
}

void goa::fitness(population &sol) {
    double tmp_best;
	for (int i = 0; i < population_num; i++){
            tmp_best=evaluation(sol[i]);
            if(tmp_best<best_obj_value){
                    best_obj_value=tmp_best;
                    best_sol=sol[i];
            }
	}
}


void goa::init(){
    // cout<<"min:"<<dmin<<endl;
    // cout<<"max:"<<dmax<<endl;
    //currentSol計得初始化;
    currentSol=population(population_num,solution(dim));
}

void goa::generate_new_solutions(){
    solution rand_centroid;
    solution tmp;
    for (int i = 0; i < population_num; i++)
        for(int j=0;j<dim;j++)
            currentSol[i][j]=frand();
    
    best_obj_value=evaluation(currentSol[0]);
    best_sol=currentSol[0];
}
double goa::evaluation(solution sol){
    double sum_of_sq=0.0;
    double sum_of_cos=0.0;
    double c=2*pi;
    for(int i=0;i<dim;i++)
        sum_of_sq +=pow(sol[i],2.0);
    for(int i=0;i<dim;i++)
        sum_of_cos +=cos(c*sol[i]);
    
    return (-20)*exp((-0.2)*sqrt((1.0/dim)*sum_of_sq))-exp((1.0/dim)*sum_of_cos)+exp(1)+20;
}
//計算兩點之間的距離
double goa::dis_cal(solution p1,solution p2){
    double sum=0.0;
    for(int i=0;i<dim;i++)
        sum+=pow(p1[i]-p2[i],2.0);
    return sum;
}

double goa::frand(){
    double f = (double)rand()/ RAND_MAX;
    return lowerbound + f * (upperbound - lowerbound);
}
double goa::coef_c(int cur_iter){
    double a;
    a=cmax-((cmax-cmin)/(double)(numIter-1))*(double)cur_iter;
    return a;
}

//----------------------------------------------------------------------force-------------------------------------------------------------------------------------//


goa::solution goa::transition(population sol,int sol_num,double c){

    double normalize_val;
    double tmp_dis;
    solution tmp_sol=solution(dim);

    for(int i=0;i<population_num;i++){
        if(sol_num==i)
            continue;
        tmp_dis=sqrt(dis_cal(sol[sol_num],sol[i]));    

        for(int j=0;j<dim;j++){
            //作正規化
            //normalize_val=((fabs(sol[i][j]-sol[sol_num][j]))/(upperbound-lowerbound))*3+1;
            if(fabs(normalize_val>2.079))
                normalize_val=fabs(normalize_val);
            else if(normalize_val==2.079)
                normalize_val=0;

            normalize_val=(tmp_dis/(upperbound-lowerbound))*3+1;
            tmp_sol[j]+=c*((upperbound-lowerbound)/2.0)*deg_forces(normalize_val)*((sol[i][j]-sol[sol_num][j])/tmp_dis);
        }    
        
    }
    //加上最好的點
    for(int i=0;i<dim;i++){

        tmp_sol[i]=c*tmp_sol[i]+best_sol[i];
        
        bool flag=true;
        while(flag){
            flag=false;
            if(tmp_sol[i]<lowerbound)
            {
                tmp_sol[i]=tmp_sol[i]+(upperbound-lowerbound)*((double)rand()/RAND_MAX);
                flag=true;
            }
            else if(tmp_sol[i]>upperbound)
            {
                tmp_sol[i]= tmp_sol[i]-(upperbound-lowerbound)*((double)rand()/RAND_MAX);
                flag=true;
            }
        } 
    }

   return tmp_sol;
}

goa::solution goa::run(){

    double tmp_best_obj_value;
    double avg_obj_value;
    population tmp_sol;
    vector<double>iter_obj_avg(numIter,0.0);
    
    init();
   
    for(int i=0;i<numRuns;i++){
        
        generate_new_solutions();
        
        fitness(currentSol);
        tmp_sol=currentSol;
        avg_obj_value=0.0;
        for(int j=0;j<numIter;j++){
            
            double c=coef_c(j);
            for(int k=0;k<population_num;k++){

                tmp_sol[k]=transition(currentSol,k,c);
                currentSol[k]=tmp_sol[k];
                /*
                tmp_best_obj_value=evaluation(tmp_sol[k]);
                
                if(tmp_best_obj_value<best_obj_value)
                {
                    best_obj_value=tmp_best_obj_value;
                    best_sol=tmp_sol[k];
                }*/
                
            }

            
             fitness(currentSol);
            //cout<<"Iteration"<<j<<":"<<best_obj_value<<endl;
            iter_obj_avg[j]+=best_obj_value;
            
        }
    }
    
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;

    cout<<"Best objective value:"<<best_obj_value<<endl;
    return best_sol;

}


#endif
       