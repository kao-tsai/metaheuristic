#ifndef __MFO_H_INCLUDED__
#define __MFO_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <random>
#include "../test_problem.h"
using namespace std;

// #define PI  3.1415926535897932384626433832795

class mfo{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		mfo(int,int,int,int,string);
		solution run();
	private:
		void init(population&);
		double evaluation(solution&);
        double frand(double,double);
        void cal_obj_and_find_best(population&);
        void update_location(int,population&,int);
        void fitness(population&);
        solution cal_mean(population&);
        
        void sort_flames(solution&,population&,int,int);
        void set_flame_no(int,population&);

private:
	int runs;
	int iters;
    int dim;
    int max_evaluate;
    double best_score;
    solution best_score_pos;
    solution cur_sol_fit;
    solution previous_sol_fit;
	population previous_sol;
    population flames;
    solution flames_fit;
    int population_num;
    solution upperbound;
    solution lowerbound;
    solution gnu_eva;
    int evaluate_num;
    string func_name;
    test_problem func;
    
    
};
mfo::mfo(int x_runs,int x_iters,int x_population_num,int x_dim,string x_func_name){
    runs=x_runs;
    iters=x_iters;
    population_num=x_population_num;
    dim=x_dim;
    
    func_name=x_func_name;
    
    func=test_problem(func_name,dim);
    
    upperbound=func.ub;
    lowerbound=func.lb;
    evaluate_num=0;
    max_evaluate=270000;
    
    
}
void mfo::init(population& current_sol){
    cur_sol_fit.assign(population_num,0);
    flames.assign(population_num,solution(dim,0));
    flames_fit.assign(population_num,0);
    evaluate_num=0;

    gnu_eva=solution(max_evaluate,0.0);
    // gnu_eva[evaluate_num]+=teacher_obj_val;
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            // cout<<j<<":"<<lowerbound[j]<<" "<<upperbound[j]<<endl;
            current_sol[i][j]=frand(lowerbound[j],upperbound[j]);
            // cout<<current_sol[i][j]<<" ";
        }
        // cout<<endl;
        cur_sol_fit[i]=evaluation(current_sol[i]);
    }
    // for(int j=0;j<dim;j++)
    //     current_sol[population_num-1][j]=0.0;
    // cout<<evaluation(current_sol[population_num-1])<<endl;
    // exit(0);

    
}

void mfo::fitness(population& all_sol){
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            if(all_sol[i][j]<lowerbound[j]){
                // all_sol[i][j]=lowerbound[j];
                // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
                all_sol[i][j]=fabs(lowerbound[j]-all_sol[i][j])+lowerbound[j];
            }
            else if(all_sol[i][j]>upperbound[j]){
                // all_sol[i][j]=upperbound[j];
                // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
                all_sol[i][j]=upperbound[j]-fabs(upperbound[j]-all_sol[i][j]);
            }
        }
    }
    for(int i=0;i<population_num;i++)
        cur_sol_fit[i]=evaluation(all_sol[i]);
    
}
double mfo::frand(double lw,double ub){
    double f = ((double)rand()/RAND_MAX);
    return lw + f * (ub-lw);
}

void mfo::update_location(int cur_iter,population &all_sol,int flame_no){
    double a=(-1.0)+(double)cur_iter*((-1.0)/iters);
    double b=1,t;
    double moth_fla_dis;
    random_device rd;
    default_random_engine gen=default_random_engine(rd());
    uniform_real_distribution<double> randn(0.0,1.0);
    for(int i=0;i<population_num;i++){
        // t=(a-1.0)*((double)rand()/RAND_MAX)+1;
        // t=(a-1.0)*randn(gen)+1;
        for(int j=0;j<dim;j++){
            if(i<=flame_no){
                moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
                // t=(a-1.0)*((double)rand()/RAND_MAX)+1;
                t=(a-1.0)*randn(gen)+1;
                all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[i][j];  
            }
            else{
                moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
                // t=(a-1.0)*((double)rand()/RAND_MAX)+1.0;
                t=(a-1.0)*randn(gen)+1;
                all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[flame_no][j];
            }
        }
    }

}
void mfo::set_flame_no(int cur_iter,population& cur_sol){
    
    if(cur_iter==1){
        // solution tmp_fit=cur_sol_fit;
        // population tmp_all_sol=cur_sol;
        // sort_flames(tmp_fit,tmp_all_sol,0,population_num-1);
        // flames=tmp_all_sol;
        // flames_fit=tmp_fit;
        sort_flames(cur_sol_fit,cur_sol,0,population_num-1);
        flames=cur_sol;
        flames_fit=cur_sol_fit;
    }
    else{
        solution double_fit=previous_sol_fit;
        population  double_all_sol=previous_sol;
        double_fit.insert(double_fit.end(),flames_fit.begin(),flames_fit.end());
        double_all_sol.insert(double_all_sol.end(),flames.begin(),flames.end());
        // cout<<double_all_sol.size()<<endl;
        sort_flames(cur_sol_fit,cur_sol,0,population_num-1);
        sort_flames(double_fit,double_all_sol,0,population_num*2-1);
        for(int i=0;i<population_num;i++){
            // cur_sol[i]=double_all_sol[i];
            // cur_sol_fit[i]=double_fit[i];
            flames[i]=double_all_sol[i];
            flames_fit[i]=double_fit[i];
        }
        
    }
    best_score=flames_fit[0];
    best_score_pos=flames[0];
    previous_sol=cur_sol;
    previous_sol_fit=cur_sol_fit;
}

void mfo::sort_flames(solution& fit,population& pop,int left,int right){
    int index;
    solution int_tmp;
    double double_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        double_tmp = fit[right];
        fit[right] = fit[index];
        fit[index] = double_tmp;

        int_tmp = pop[right];
        pop[right] = pop[index];
        pop[index] = int_tmp;

        pivot = fit[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (fit[j] <= pivot)
            {
                i+=1;
                double_tmp = fit[j];
                fit[j] = fit[i];
                fit[i] = double_tmp;
                
                int_tmp = pop[j];
                pop[j] = pop[i];
                pop[i] = int_tmp;
            }
        }
        index=i+1;
        double_tmp = fit[index];
        fit[index] = fit[right];
        fit[right] = double_tmp;

        int_tmp = pop[index];
        pop[index] = pop[right];
        pop[right] = int_tmp;

        sort_flames(fit,pop, left, index-1);
        sort_flames(fit,pop, index+1, right);
    }
}

mfo::solution mfo::run(){
    // solution try_sol(dim,0.0);
    // cout<<evaluation(try_sol)<<endl; 
    // exit(0);
    population current_sol(population_num,solution(dim));
    solution each_run_best;
    double each_run_best_obj=numeric_limits<double>::infinity();;
    int flame_no;
    vector<double>iter_obj_avg(iters,0.0);
    
    for(int i=0;i<runs;i++){
        
        init(current_sol);
        // for(int k=0;k<population_num;k++){
        //     cout<<k<<" "<<current_sol[k][0]<<"--->"<<cur_sol_fit[k]<<endl;
        // }
        // sort_flames(cur_sol_fit,current_sol,0,population_num-1);
        // for(int k=0;k<population_num;k++){
        //     cout<<k<<" "<<current_sol[k][0]<<"--->"<<cur_sol_fit[k]<<endl;
        // }
        // exit(0);
        
        for(int j=0;j<iters;j++){
            
            flame_no=round(population_num-j*((double)(population_num-1)/(double)iters));
            
            fitness(current_sol);

            set_flame_no(j+1,current_sol);

            update_location(j+1,current_sol,flame_no);
            
            
            iter_obj_avg[j]+=best_score;
        }
        if(best_score<each_run_best_obj){
            each_run_best_obj=best_score;
            each_run_best=best_score_pos;
        }
        
    }
    ofstream output("result/mfo/"+func_name+".txt");
    for(int i=0;i<iters;i++)
        output<<i<<" "<<(double)iter_obj_avg[i]/runs<<endl;
    output.close();
    for(int i=0;i<iters;i++)
        cout<<fixed<<setprecision(20)<<(double)iter_obj_avg[i]/runs<<endl;
    // cout<<"Evaluation Number:"<<evaluate_num<<endl;
    cout<<"Best objective value:"<<each_run_best_obj<<endl;
    return  each_run_best;

}

double mfo::evaluation(solution& sol){
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