#ifndef __HOO_H_INCLUDED__
#define __HOO_H_INCLUDED__
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
class hoo{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		hoo(int,int,int);
		solution run();
	private:
		void init(population&);
		double evaluation(solution&);
        double frand();
        void update_location(population&,double);
        void determination(population&);
        solution cal_mean(population&);
        solution levy_flight();

private:
	int runs;
	int iters;
    int dim;
    int max_evaluate;
	double rabbit_obj_val;//Best objective value
	solution rabbit_location;//Best location of rabbit
    int population_num;
    double upperbound;
    double lowerbound;
    solution gnu_eva;
    int evaluate_num;
};
hoo::hoo(int x_runs,int x_iters,int x_population_num){
    runs=x_runs;
    iters=x_iters;
    population_num=x_population_num;
    upperbound=30;
    lowerbound=-30;
    dim=30;
    evaluate_num=0;
    max_evaluate=270000;
}
void hoo::determination(population& all_sol){
    double fitness;
    //檢查每一組解是否超過邊界 以及 更新兔子的位置和最佳值
    for(int i=0;i<population_num;i++){
        for(int j=0;j<dim;j++){
            if(all_sol[i][j]>upperbound)
                all_sol[i][j]=upperbound;
            else if(all_sol[i][j]<lowerbound)
                all_sol[i][j]=lowerbound;
        }
        fitness=evaluation(all_sol[i]);
        if(fitness<rabbit_obj_val){
            rabbit_obj_val=fitness;
            rabbit_location=all_sol[i];
        }
    }
}
double hoo::frand(){
    double f = (double)rand()/ RAND_MAX;
    return lowerbound + f * (upperbound - lowerbound);
}

double hoo::evaluation(solution& sol){
    double sum_of_sq=0.0;
    double sum_of_cos=0.0;
    double c=2*pi;
    for(int i=0;i<dim;i++)
        sum_of_sq +=pow(sol[i],2.0);
    for(int i=0;i<dim;i++)
        sum_of_cos +=cos(c*sol[i]);
    
    return (-20)*exp((-0.2)*sqrt((1.0/dim)*sum_of_sq))-exp((1.0/dim)*sum_of_cos)+exp(1)+20;
}


// double hoo::evaluation(solution& sol){
//     double sum=0.0;
//     for(int i=0;i<dim-1;i++)
//         sum+=100.0*pow(sol[i+1]-pow(sol[i],2.0),2.0)+pow(sol[i]-1,2.0);
//     evaluate_num++;
//     gnu_eva[evaluate_num]+=rabbit_obj_val;
//     return sum;
// }

void hoo::update_location(population &hawks,double E1){
    double E0,escaping_energy,q;
    for(int i=0;i<population_num;i++){
        E0=2.0*((double)rand()/RAND_MAX)-1;
        escaping_energy=E0*E1;
        if(fabs(escaping_energy)>=1){
            q=(double)rand()/RAND_MAX;
            if(q<0.5){
                solution rand_hawk=hawks[rand()%population_num];
                double r1=(double)rand()/RAND_MAX;
                double r2=(double)rand()/RAND_MAX;
                for(int j=0;j<dim;j++)
                    hawks[i][j]=rand_hawk[j]-r1*fabs(rand_hawk[j]-2.0*r2*hawks[i][j]); 
            }
            else if(q>=0.5){
                double r1=(double)rand()/RAND_MAX;
                double r2=(double)rand()/RAND_MAX;
                solution mean=cal_mean(hawks);
                
                for(int j=0;j<dim;j++)
                    hawks[i][j]=(rabbit_location[j]-mean[j])-r1*((upperbound-lowerbound)*r2+lowerbound);
            }
        }
        else if(fabs(escaping_energy)<1){
            //do the exploitation phase
            double r=(double)rand()/RAND_MAX;
            if(r>=0.5 && fabs(escaping_energy)<0.5){
                for(int j=0;j<dim;j++)
                    hawks[i][j]=rabbit_location[j]-escaping_energy*fabs(rabbit_location[j]-hawks[i][j]);
            }
            if(r>=0.5 && fabs(escaping_energy)>=0.5){
                double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                for(int j=0;j<dim;j++)
                    hawks[i][j]=(rabbit_location[j]-hawks[i][j])-escaping_energy*fabs(jump_strength*rabbit_location[j]-hawks[i][j]);
            }
            if(r<0.5 && fabs(escaping_energy)>=0.5){
                //
                double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                solution get_levy=levy_flight();
                solution X1(dim),X2(dim);
                for(int j=0;j<dim;j++)
                    X1[j]=rabbit_location[j]-escaping_energy*fabs(jump_strength*rabbit_location[j]-hawks[i][j]);
                if(evaluation(X1)<evaluation(hawks[i]))
                    hawks[i]=X1;
                else
                {
                    for(int j=0;j<dim;j++)
                        X2[j]=X1[j]+((double)rand()/RAND_MAX)*get_levy[j];
                    if(evaluation(X2)<evaluation(hawks[i]))
                        hawks[i]=X2;
                }

            }
            if(r<0.5 && fabs(escaping_energy)<0.5){
                //
                double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                solution get_levy=levy_flight();
                solution X1(dim),X2(dim);
                solution mean=cal_mean(hawks);
                for(int j=0;j<dim;j++)
                    X1[j]=rabbit_location[j]-escaping_energy*fabs(jump_strength*rabbit_location[j]-mean[j]);
                if(evaluation(X1)<evaluation(hawks[i]))
                    hawks[i]=X1;
                else
                {
                    for(int j=0;j<dim;j++)
                        X2[j]=X1[j]+((double)rand()/RAND_MAX)*get_levy[j];//Debug:多乘以1個levy flight;
                    if(evaluation(X2)<evaluation(hawks[i]))
                        hawks[i]=X2;
                }
            }
        }
    }
}
hoo::solution hoo::levy_flight(){
    random_device rd;
    default_random_engine gen=default_random_engine(rd());
    uniform_real_distribution<double> randn(0.0,1.0);
    // normal_distribution<double> randn(0.0,1.0);
    double beta=1.5;
    double sigma=(gamma(1+beta)*sin((beta/2.0)*M_PI)/(gamma((1+beta)/2)*beta*pow((beta-1)/2,2)));
    sigma=pow(sigma,(1/beta));
    solution u(dim);
    for(int i=0;i<dim;i++)
        u[i]=(randn(gen)*sigma)/pow(fabs(randn(gen)),1.0/beta);
    return u;
}

hoo::solution hoo::cal_mean(population& all_sol){
    solution sum(dim,0.0);
    for(int i=0;i<population_num;i++)
        transform(sum.begin(),sum.end(),all_sol[i].begin(),sum.begin(),plus<double>());
    for(int i=0;i<dim;i++)
        sum[i]/=population_num;
    return sum;
}

void hoo::init(population& current_sol){

    evaluate_num=0;
    rabbit_obj_val=numeric_limits<double>::infinity();
    gnu_eva=solution(max_evaluate,0.0);
    gnu_eva[evaluate_num]+=rabbit_obj_val;
    for(int i=0;i<population_num;i++)
        for(int j=0;j<dim;j++)
            current_sol[i][j]=frand();
    rabbit_location=solution(dim,0.0);
    
}

hoo::solution hoo::run(){
    population current_sol(population_num,solution(dim));
    vector<double>iter_obj_avg(iters,0.0);
    double E1;
    for(int i=0;i<runs;i++){
        init(current_sol);
        for(int j=0;j<iters;j++){
            determination(current_sol);
            E1=2.0*(1.0-((double)j/(double)iters));
            update_location(current_sol,E1);
            iter_obj_avg[j]+=rabbit_obj_val;
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
    cout<<"Best objective value:"<<rabbit_obj_val<<endl;
    return  rabbit_location;

}



#endif