#ifndef __PSO_H_INCLUDED__
#define __PSO_H_INCLUDED__
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
class pso{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		pso(int,int,int,int,string);
		solution run();
	private:
		void init(population&);
        double frand(double , double);
        void update_location(population&,population&);
        void determination(population&);
        double evaluation(solution&);
        // solution cal_mean(population&);
        // solution levy_flight();

private:
	int runs;
	int iters;
    int dim;
    int max_evaluate;
	double gbest_obj;//Best objective value
	solution gbest_pos;//Best location of rabbit
    population pbest_pos;
    solution pbest_obj;
    int population_num;
    solution upperbound;
    solution lowerbound;
    solution gnu_eva;
    int evaluate_num;
    test_problem func;
    string func_name;
};
pso::pso(int x_runs,int x_iters,int x_population_num,int x_dim,string x_func_name){
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
void pso::determination(population& all_sol){
    double fitness;

    //檢查每一組解是否超過邊界 以及 更新兔子的位置和最佳值
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            if(all_sol[i][j]>upperbound[j]){
                // all_sol[i][j]=upperbound[j];
                all_sol[i][j]=fabs(lowerbound[j]-all_sol[i][j])+lowerbound[j];
            }
            else if(all_sol[i][j]<lowerbound[j]){
                all_sol[i][j]=upperbound[j]-fabs(upperbound[j]-all_sol[i][j]);
                // all_sol[i][j]=lowerbound[j];
            }
        }
        fitness=evaluation(all_sol[i]);
        if(fitness<pbest_obj[i]){
            pbest_obj[i]=fitness;
            pbest_pos[i]=all_sol[i];
        }
        if(fitness<gbest_obj){
            gbest_obj=fitness;
            gbest_pos=all_sol[i];
        }
    }
}
double pso::frand(double lb,double ub){
    double f = (double)rand()/ RAND_MAX;
    return lb + f * (ub - lb);
}

void pso::update_location(population &all_sol,population& Vij){
    double w=0.5,c1=2,c2=2;
    double r1,r2;
    
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            r1=(double)rand()/RAND_MAX;
            r2=(double)rand()/RAND_MAX;
            Vij[i][j]=w*Vij[i][j]+
            c1*r1*(pbest_pos[i][j]-all_sol[i][j])+c2*r2*(gbest_pos[j]-all_sol[i][j]);
            all_sol[i][j]=all_sol[i][j]+Vij[i][j];
        }
    }

}

void pso::init(population& current_sol){

    evaluate_num=0;
    gbest_obj=numeric_limits<double>::infinity();
    pbest_pos.assign(population_num,solution(dim,0));
    pbest_obj.assign(population_num,0);
    gnu_eva=solution(max_evaluate,0.0);
    gnu_eva[evaluate_num]+=gbest_obj;
    for(int i=0;i<population_num;i++)
        for(int j=0;j<dim;j++)
            current_sol[i][j]=frand(lowerbound[j],upperbound[j]);
    gbest_pos=solution(dim,0.0);
    for(int i=0;i<population_num;i++){
        pbest_pos[i]=current_sol[i];
        pbest_obj[i]=evaluation(current_sol[i]);
        if(pbest_obj[i]<gbest_obj){
            gbest_obj=pbest_obj[i];
            gbest_pos=pbest_pos[i];
        }
    }
    
}

pso::solution pso::run(){
    population current_sol(population_num,solution(dim));
    solution each_run_best;
    double each_run_best_obj=numeric_limits<double>::infinity();;
    vector<double>iter_obj_avg(iters,0.0);
    population vel;
    for(int i=0;i<runs;i++){
        init(current_sol);
        vel=population(population_num,solution(dim,0));
        for(int j=0;j<iters;j++){
            update_location(current_sol,vel);
            determination(current_sol);
            iter_obj_avg[j]+=gbest_obj;
            
        }
        if(gbest_obj<each_run_best_obj){
            each_run_best_obj=gbest_obj;
            each_run_best=gbest_pos;
        }
    }
    
    ofstream output("result/pso/"+func_name+".txt");
    for(int i=0;i<iters;i++)
        output<<i<<" "<<(double)iter_obj_avg[i]/runs<<endl;
    output.close();
    for(int i=0;i<iters;i++)
        cout<<fixed<<setprecision(20)<<(double)iter_obj_avg[i]/runs<<endl;
    // cout<<"Evaluation Number:"<<evaluate_num<<endl;
    cout<<"Best objective value:"<<each_run_best_obj<<endl;
    return  each_run_best;

}

double pso::evaluation(solution& sol){
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