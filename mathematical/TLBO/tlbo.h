#ifndef __TLBO_H_INCLUDED__
#define __TLBO_H_INCLUDED__
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
#define pi 3.141592653;

class tlbo{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		tlbo(int,int,int,int,string);
		solution run();
	private:
		void init(population&);
		double evaluation(solution&);
        double frand(double,double);
        void cal_obj_and_find_best(population&);
        void update_location(population&);
        void determination(population&);
        solution cal_mean(population&);
        // solution levy_flight();
        void find_teacher(population&);

private:
	int runs;
	int iters;
    int dim;
    int max_evaluate;
    solution all_sol_obj;
	double teacher_obj_val;//Best objective value
	solution teacher_location;//Best location of rabbit
    int population_num;
    solution upperbound;
    solution lowerbound;
    solution gnu_eva;
    int evaluate_num;
    string func_name;
    test_problem func;
};
tlbo::tlbo(int x_runs,int x_iters,int x_population_num,int x_dim,string x_func_name){
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

void tlbo::determination(population& all_sol){
    double temp_best_obj=numeric_limits<double>::infinity();
    solution temp_best_sol;
    for(int i=0;i<population_num;i++){
        if(temp_best_obj>all_sol_obj[i]){
            temp_best_obj=all_sol_obj[i];
            temp_best_sol=all_sol[i];
        }
    }
    teacher_location=temp_best_sol;
    teacher_obj_val=temp_best_obj;
}
double tlbo::frand(double lw,double ub){
    double f = (double)rand()/ RAND_MAX;
    return lw + f * (ub-lw);
}
/*
double tlbo::evaluation(solution& sol){
    double sum_of_sq=0.0;
    double sum_of_cos=0.0;
    double c=2*pi;
    for(int i=0;i<dim;i++)
        sum_of_sq +=pow(sol[i],2.0);
    for(int i=0;i<dim;i++)
        sum_of_cos +=cos(c*sol[i]);
    
    return (-20)*exp((-0.2)*sqrt((1.0/dim)*sum_of_sq))-exp((1.0/dim)*sum_of_cos)+exp(1)+20;
}
*/
/*
double tlbo::evaluation(solution& sol){
    double sum=0.0;
    for(int i=0;i<dim-1;i++)
        sum+=100.0*pow(sol[i+1]-pow(sol[i],2.0),2.0)+pow(sol[i]-1,2.0);
    evaluate_num++;
    gnu_eva[evaluate_num]+=teacher_obj_val;
    return sum;
}*/

void tlbo::update_location(population &all_sol){
    solution old_mean=cal_mean(all_sol);
    double TF=rand()%2+1;
    population new_sol(population_num,solution(dim));
    solution new_sol_obj(population_num);
    int Sj;
    //----------------------------Teacher phase-----------------------------//
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            //Get new location
            new_sol[i][j]=all_sol[i][j]+((double)rand()/RAND_MAX)*(teacher_location[j]-(TF)*old_mean[j]);
            
            //Keep the value in the between the upperbound and lowerbound
            if(new_sol[i][j]>upperbound[j])
                new_sol[i][j]=upperbound[j];
            else if(new_sol[i][j]<lowerbound[j])
                new_sol[i][j]=lowerbound[j];
        }
        new_sol_obj[i]=evaluation(new_sol[i]);
    
    
        //replace solution if new solution is better than old solution  
        if(new_sol_obj[i]<all_sol_obj[i]){
            all_sol[i]=new_sol[i];
            all_sol_obj[i]=new_sol_obj[i];
        }
    
    //-------------------------------------------------------------------------//
    //-----------------------------Student Phase-------------------------------//
    
    // for(int i=0;i<population_num;i++){
        do{
            Sj=rand()%population_num;
        }while(Sj==i);
        for(int j=0;j<dim;j++){
            if(all_sol_obj[i]<all_sol_obj[Sj])
                new_sol[i][j]=all_sol[i][j]+((double)rand()/RAND_MAX)*(all_sol[i][j]-all_sol[Sj][j]);
            else
                new_sol[i][j]=all_sol[i][j]+((double)rand()/RAND_MAX)*(all_sol[Sj][j]-all_sol[i][j]);
            //Keep the value in the between the upperbound and lowerbound
            if(new_sol[i][j]>upperbound[j])
                new_sol[i][j]=upperbound[j];
            else if(new_sol[i][j]<lowerbound[j])
                new_sol[i][j]=lowerbound[j];
        }
    // }
    //replace solution if new solution is better than old solution 
    if(new_sol_obj[i]<all_sol_obj[i]){
        all_sol[i]=new_sol[i];
        all_sol_obj[i]=new_sol_obj[i];
    }
    //-------------------------------------------------------------------------//
    }

}
//calculate each solution's objective value and find the best solution as teacher
void tlbo::cal_obj_and_find_best(population& all_sol){
    double temp_best_obj=numeric_limits<double>::infinity();
    solution temp_best_sol;
    for(int i=0;i<population_num;i++){
        all_sol_obj[i]=evaluation(all_sol[i]);
        if(temp_best_obj>all_sol_obj[i]){
            temp_best_obj=all_sol_obj[i];
            temp_best_sol=all_sol[i];
        }
    }
    teacher_location=temp_best_sol;
    teacher_obj_val=temp_best_obj;
}
/*
tlbo::solution tlbo::levy_flight(){
    default_random_engine gen;
    normal_distribution<double> randn(0.0,1.0);
    double beta=1.5;
    double sigma=(gamma(1+beta)*sin((beta/2.0)*M_PI)/(gamma((1+beta)/2)*beta*pow((beta-1)/2,2)));
    sigma=pow(sigma,(1/beta));
    solution u(dim);
    for(int i=0;i<dim;i++)
        u[i]=(randn(gen)*sigma)/pow(fabs(randn(gen)),1.0/beta);
    return u;
}
*/

tlbo::solution tlbo::cal_mean(population& all_sol){
    solution sum(dim,0.0);
    for(int i=0;i<population_num;i++)
        transform(sum.begin(),sum.end(),all_sol[i].begin(),sum.begin(),plus<double>());
    for(int i=0;i<dim;i++)
        sum[i]/=population_num;
    return sum;
}

void tlbo::init(population& current_sol){
    all_sol_obj.assign(population_num,-1);
    evaluate_num=0;
    teacher_obj_val=numeric_limits<double>::infinity();
    gnu_eva=solution(max_evaluate,0.0);
    gnu_eva[evaluate_num]+=teacher_obj_val;
    for(int i=0;i<population_num;i++)
        for(int j=0;j<dim;j++)
            current_sol[i][j]=frand(lowerbound[j],upperbound[j]);
    teacher_location=solution(dim,-1);
    
}

tlbo::solution tlbo::run(){
    population current_sol(population_num,solution(dim));
    vector<double>iter_obj_avg(iters,0.0);
    
    for(int i=0;i<runs;i++){
        
        init(current_sol);
        
        cal_obj_and_find_best(current_sol);
        
        for(int j=0;j<iters;j++){
            //Modify solutions based on best solution
            update_location(current_sol);
            
            determination(current_sol);
            iter_obj_avg[j]+=teacher_obj_val;
        }
    }
    ofstream output("gnuplot/result/HHO.txt");
    for(int i=0;i<max_evaluate;i++)
        if(i%500==0){
            output<<i<<" ";
            output<<fixed<<setprecision(20)<<gnu_eva[i]/runs<<endl;
        }
    output.close();
    for(int i=0;i<iters;i++)
        cout<<fixed<<setprecision(20)<<(double)iter_obj_avg[i]/runs<<endl;
    cout<<"Evaluation Number:"<<evaluate_num<<endl;
    cout<<"Best objective value:"<<teacher_obj_val<<endl;
    return  teacher_location;

}

double tlbo::evaluation(solution& sol){
    if(func_name=="Ackley")
        return func.Ackley(sol);
    
}

#endif