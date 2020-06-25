#ifndef __DE_H_INCLUDED__
#define __DE_H_INCLUDED__
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

class de{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		de(int,int,int,int,string);
		solution run();
	private:
		void init(population&);
		double evaluation(solution&);
        double frand(double,double);
        void cal_obj_and_find_best(population&);
        void update_location(population&);
        void fitness(population&);
        solution cal_mean(population&);
        void quick_sort(solution&,vector<int>&,int ,int);

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
    double mutation_pro;
    solution gbest;
    double gbest_fit;
    double inf=numeric_limits<double>::infinity();
    
};
de::de(int x_runs,int x_iters,int x_population_num,int x_dim,string x_func_name){
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
    mutation_pro=1.0/dim;
    
}
void de::init(population& current_sol){
    cur_sol_fit.assign(population_num,0);
    flames.assign(population_num,solution(dim,0));
    flames_fit.assign(population_num,0);
    evaluate_num=0;
    gbest_fit=inf;
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
        if(cur_sol_fit[i]<gbest_fit){
            gbest_fit=cur_sol_fit[i];
            gbest=current_sol[i];
        }
            
    }
    
    
    
}

void de::fitness(population& all_sol){
    // for(int i=0;i<population_num;i++){
    //     for(int j=0;j<dim;j++){
    //         if(all_sol[i][j]<lowerbound[j]){
    //             // all_sol[i][j]=lowerbound[j];
    //             // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
    //             all_sol[i][j]=fabs(lowerbound[j]-all_sol[i][j])+lowerbound[j];
    //         }
    //         else if(all_sol[i][j]>upperbound[j]){
    //             // all_sol[i][j]=upperbound[j];
    //             // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
    //             all_sol[i][j]=upperbound[j]-fabs(upperbound[j]-all_sol[i][j]);
    //         }
    //     }
    // }
    for(int i=0;i<population_num;i++)
        cur_sol_fit[i]=evaluation(all_sol[i]);
    
}
double de::frand(double lw,double ub){
    double f = ((double)rand()/RAND_MAX);
    return lw + f * (ub-lw);
}

void de::update_location(population &all_sol){
    
    double F=0.2;
    double cr_rate=0.7;
    population new_sol=all_sol;
    random_device rd;
    default_random_engine gen=default_random_engine(rd());
    uniform_real_distribution<double> randn(0.0,1.0);
    for(int i=0;i<population_num;i++){
       

        double r2=rand()%population_num;
        double r3=rand()%population_num;
        do{
            r2=rand()%population_num;
        }while(r2==i);

        do{
            r2=rand()%population_num;
        }while(r2==i&&r3==r2);

        for(int j=0;j<dim;j++){
            double rnd=(double)rand()/RAND_MAX;
 
            if(rnd<cr_rate){
                new_sol[i][j]=all_sol[i][j]+F*(gbest[j]-new_sol[i][j])+F*(all_sol[r2][j]-all_sol[r3][j]);
            }
            else{
               new_sol[i][j]=all_sol[i][j];
            }
            if(new_sol[i][j]<lowerbound[j]){
            // all_sol[i][j]=lowerbound[j];
            // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
            new_sol[i][j]=fabs(lowerbound[j]-all_sol[i][j])+lowerbound[j];
            }
            else if(new_sol[i][j]>upperbound[j]){
                // all_sol[i][j]=upperbound[j];
                // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
                new_sol[i][j]=upperbound[j]-fabs(upperbound[j]-all_sol[i][j]);
            }
            // if(isnan(all_sol[i][j]))
            //     all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
            double new_fit=evaluation(new_sol[i]);
            if(new_fit<cur_sol_fit[i]){
                all_sol[i]=new_sol[i];
                cur_sol_fit[i]=new_fit;
            }
            if(new_fit<gbest_fit){
                gbest=new_sol[i];
                gbest_fit=new_fit;
            }

        }
    }
    

}


void de::quick_sort(solution& value,vector<int>& this_rank,int left,int right){
    int index;
    int int_tmp;
    double double_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        double_tmp = value[right];
        value[right] = value[index];
        value[index] = double_tmp;

        int_tmp = this_rank[right];
        this_rank[right] = this_rank[index];
        this_rank[index] = int_tmp;

        pivot = value[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (value[j] <= pivot)
            {
                i+=1;
                double_tmp = value[j];
                value[j] = value[i];
                value[i] = double_tmp;
                
                int_tmp = this_rank[j];
                this_rank[j] = this_rank[i];
                this_rank[i] = int_tmp;
            }
        }
        index=i+1;
        double_tmp = value[index];
        value[index] = value[right];
        value[right] = double_tmp;

        int_tmp = this_rank[index];
        this_rank[index] = this_rank[right];
        this_rank[right] = int_tmp;

        quick_sort(value,this_rank, left, index-1);
        quick_sort(value,this_rank, index+1, right);
    }
}

de::solution de::run(){
    
    population current_sol(population_num,solution(dim));
    solution each_run_best;
    double each_run_best_obj=numeric_limits<double>::infinity();;
    int flame_no;
    double q_r;
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
            
            fitness(current_sol);   
            
            
           
            update_location(current_sol);
            
            

            // for(int k=0;k<population_num;k++){
            //     if(cur_sol_fit[k]>=previous_sol_fit[k]){
                    
            //     }
            // }
            
            
            
            iter_obj_avg[j]+=best_score;
        }
        if(gbest_fit<each_run_best_obj){
            each_run_best_obj=gbest_fit;
            each_run_best=gbest;
        }
        // cout<<"run:"<<i<<" "<<best_score<<endl;
    }
    ofstream output("result/de/"+func_name+".txt");
    for(int i=0;i<iters;i++)
        output<<i<<" "<<(double)iter_obj_avg[i]/runs<<endl;
    output.close();
    // for(int i=0;i<iters;i++)
    //     cout<<(double)iter_obj_avg[i]/runs<<endl;

    cout<<func_name<<": "<<(double)iter_obj_avg[iters-1]/runs<<endl;

    // cout<<"Evaluation Number:"<<evaluate_num<<endl;
    // cout<<"Best objective value:"<<each_run_best_obj<<endl;
    return  each_run_best;

}

double de::evaluation(solution& sol){
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
    else if(func_name=="Griewank10")
        return func.Griewank10(sol);
    else if(func_name=="Schaffer2")
        return func.Schaffer2(sol);
    else if(func_name=="Schwefel")
        return func.Schwefel(sol);
    else if(func_name=="Bohachevsky1")
        return func.Bohachevsky1(sol);
    else if(func_name=="Sum_Square")
        return func.Sum_Square(sol);
    else if(func_name=="Booth")
        return func.Booth(sol);
    else if(func_name=="Zakharov")
        return func.Zakharov(sol);
    else if(func_name=="Three_Hump_Camel")
        return func.Three_Hump_Camel(sol);
    else if(func_name=="De_Jong5")
        return func.De_Jong5(sol);
    else if(func_name=="Beale")
        return func.Beale(sol);
    else if(func_name=="Powell")
        return func.Powell(sol);

}

#endif