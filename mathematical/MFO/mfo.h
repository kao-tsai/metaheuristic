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
        population clonal_selection();
        void sort_flames(solution&,population&,int,int);
        void set_flame_no(int,population&);
        void mutation(solution&);
        void clonal_selection(population&);

        void generate_new_Cr_F();
        void update_uCr_uF();
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
    double u_b;
    double u_t;
    solution set_b;
    solution set_t;
    solution self_b;
    solution self_t;
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
    mutation_pro=1.0/dim;
    
}
void mfo::init(population& current_sol){
    cur_sol_fit.assign(population_num,0);
    flames.assign(population_num,solution(dim,0));
    flames_fit.assign(population_num,0);
    evaluate_num=0;
    u_b=0.2;
    u_t=0.2;
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
    self_b.assign(population_num,0.2);
    self_t.assign(population_num,0.2);
    set_b.clear();
    set_t.clear();
    
}

void mfo::fitness(population& all_sol){
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
double mfo::frand(double lw,double ub){
    double f = ((double)rand()/RAND_MAX);
    return lw + f * (ub-lw);
}

void mfo::update_location(int cur_iter,population &all_sol,int flame_no){
    double a=(-1.0)+(double)cur_iter*((-1.0)/iters);
    double b=1,t;
    double moth_fla_dis;
    double F=0.2;
    double cr_rate=0.5;
    random_device rd;
    default_random_engine gen=default_random_engine(rd());
    uniform_real_distribution<double> randn(0.0,1.0);
    for(int i=0;i<population_num;i++){
        // t=(a-1.0)*((double)rand()/RAND_MAX)+1;
        // t=(a-1.0)*randn(gen)+1;
        // for(int j=0;j<dim;j++){
        //     if(i<=flame_no){
        //         moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
        //         // t=(a-1.0)*((double)rand()/RAND_MAX)+1;
        //         t=(a-1.0)*randn(gen)+1;
        //         all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[i][j];
        //         // all_sol[i][j]=moth_fla_dis*exp(self_b[i]*t)*cos(t*2*PI)+flames[i][j];
               
        //     }
        //     else{
        //         moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
        //         // t=(a-1.0)*((double)rand()/RAND_MAX)+1.0;
        //         t=(a-1.0)*randn(gen)+1;
        //         all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[flame_no][j];
        //         // all_sol[i][j]=moth_fla_dis*exp(self_b[i]*t)*cos(t*2*PI)+flames[flame_no][j];
                
        //     }
        // }

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
            if(i<=flame_no){
                if(rnd<cr_rate){
                    all_sol[i][j]=all_sol[i][j]+F*(flames[i][j]-all_sol[i][j])+F*(all_sol[r2][j]-all_sol[r3][j]);
                }
                else{
                    moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
                    // t=(a-1.0)*((double)rand()/RAND_MAX)+1;
                    t=(a-1.0)*randn(gen)+1;
                    all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[i][j];
                }
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
                if(isnan(all_sol[i][j]))
                    all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
            }
            else{
                if(rnd<cr_rate){
                    all_sol[i][j]=all_sol[i][j]+F*(flames[flame_no][j]-all_sol[i][j])+F*(all_sol[r2][j]-all_sol[r3][j]);
                }
                else{
                moth_fla_dis=fabs(flames[i][j]-all_sol[i][j]);
                // t=(a-1.0)*((double)rand()/RAND_MAX)+1.0;
                t=(a-1.0)*randn(gen)+1;
                all_sol[i][j]=moth_fla_dis*exp(b*t)*cos(t*2*PI)+flames[flame_no][j];
                }
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
    }
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

}

void mfo::generate_new_Cr_F(){
    random_device rd;
    mt19937 generator(rd());
    normal_distribution<double> b_dist(u_b,0.1);
    normal_distribution<double> t_dist(u_t,0.1);

    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1);
    
    int one_third=population_num/3;
    for(int i=0;i<population_num;i++)
    {

            self_b[i]=b_dist(generator);
            self_t[i]=t_dist(generator);
            // if(((double)rand()/RAND_MAX)<0.1)
            //     self_b[i]=0.0+((double)rand()/RAND_MAX)*(0.9-0.0);
            
            // if(((double)rand()/RAND_MAX)<0.1)
            //     self_t[i]=0.1+((double)rand()/RAND_MAX)*(0.9-0.1);
    }
    // for(int i=0;i<one_third;i++){
    //     int choose=rand()%population_num;
    //     self_t[choose]=rndF(gen);
    // }
    set_b.clear();
    set_t.clear();
}

void mfo::update_uCr_uF(){
    double c=0.1;
    //
    double mean_S_cr=0;
    for(int i=0;i<set_b.size();i++)
        mean_S_cr+=set_b[i];
    if(mean_S_cr!=0)
        mean_S_cr/=set_b.size();
    u_b=(1-c)*u_b+c*mean_S_cr;
    //
    double lemer=0.0;
    double sum1=0,sum2=0;
    for(int i=0;i<set_t.size();i++){
        sum1+=pow(set_t[i],2.0);
        sum2+=set_t[i];
    }
    if(sum2!=0)
        lemer=sum1/sum2;
    u_t=(1-c)*u_t+c*lemer;
}
void mfo::clonal_selection(population& all_sol){
    int NC=6;
    population PC(NC);
    solution PC_fit(NC);

    solution double_fit=previous_sol_fit;
    population  double_all_sol=previous_sol;
    double_fit.insert(double_fit.end(),flames_fit.begin(),flames_fit.end());
    double_all_sol.insert(double_all_sol.end(),flames.begin(),flames.end());
    // cout<<double_all_sol.size()<<endl;
    sort_flames(cur_sol_fit,all_sol,0,population_num-1);
    sort_flames(double_fit,double_all_sol,0,population_num*2-1);
    for(int i=0;i<NC;i++){
        PC[i]=double_all_sol[i];
        PC_fit[i]=double_fit[i];
    }

    // for(int i=0;i<population_num;i++){
    //     // cur_sol[i]=double_all_sol[i];
    //     // cur_sol_fit[i]=double_fit[i];
    //     flames[i]=double_all_sol[i];
    //     flames_fit[i]=double_fit[i];
    // }

    population new_pop_E;
    solution new_pop_fit;
    double fit_sum=0.0;
    for(int i=0;i<NC;i++)
        fit_sum+=PC_fit[i];

    int q_i;
    bool isfull=false;
    int test_sum=0;
    for(int i=0;i<NC;i++){
        q_i=ceil(((double)population_num)*(PC_fit[i]/fit_sum));
       
        for(int j=0;j<q_i;j++){
            new_pop_E.push_back(PC[i]);
            new_pop_fit.push_back(PC_fit[i]);
            if(new_pop_E.size()==population_num){
                isfull=true;
                break;
            }
        }
        if(isfull)
            break;
    }
    flames=new_pop_E;
    flames_fit=new_pop_fit;
    best_score=flames_fit[0];
    best_score_pos=flames[0];
    previous_sol=all_sol;
    previous_sol_fit=cur_sol_fit;

}

void mfo::mutation(solution& x1){
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    double dimm=20.0;
    random_device rd;
    mt19937 generator(rd());
    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1);
    for (j=0; j<dim; j++)
    {
        rnd=rndF(gen);
        if (rnd <= mutation_pro)
        {   
           
            y =x1[j];

            yl = lowerbound[j];
            yu = upperbound[j];
            delta1 = (y-yl)/(yu-yl);

            delta2 = (yu-y)/(yu-yl);
            rnd = rndF(gen);
            mut_pow = 1.0/(dimm+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(dimm+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(dimm+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            x1[j] = y;

            if(x1[j]<lowerbound[j]){
                // all_sol[i][j]=lowerbound[j];
                // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
                x1[j]=fabs(lowerbound[j]-x1[j])+lowerbound[j];
            }
            else if(x1[j]>upperbound[j]){
                // all_sol[i][j]=upperbound[j];
                // all_sol[i][j]=frand(lowerbound[j],upperbound[j]);
                x1[j]=upperbound[j]-fabs(upperbound[j]-x1[j]);
            }
        }
    }

}

void mfo::set_flame_no(int cur_iter,population& cur_sol){
    vector<int>ori_id(population_num);
    for(int i=0;i<population_num;i++)
        ori_id[i]=i;
    
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
        
        // quick_sort(cur_sol_fit,ori_id,0,population_num-1);
        // solution sort_b=self_b;
        // solution sort_t=self_t;
        // population sort_pop=cur_sol;
        // for(int i=0;i<population_num;i++)
        // {
        //     self_b[i]=sort_b[ori_id[i]];
        //     self_t[i]=self_t[ori_id[i]];
        //     cur_sol[i]=sort_pop[ori_id[i]];
        // }
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
void mfo::quick_sort(solution& value,vector<int>& this_rank,int left,int right){
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
            
            flame_no=round(population_num-j*((double)(population_num-1)/(double)iters));
            fitness(current_sol);
            // if(j!=0){
            // for(int k=0;k<population_num;k++){
            //     if(cur_sol_fit[k]<previous_sol_fit[k]){
            //         set_b.push_back(self_b[k]);
            //         set_t.push_back(self_t[k]);
            //     }
            // }
            // }
            // // update_uCr_uF();

            
            // q_r=(double)rand()/RAND_MAX;
            // if(j<80000)
                set_flame_no(j+1,current_sol);
            // else
            //     clonal_selection(current_sol);
            // generate_new_Cr_F();
            update_location(j+1,current_sol,flame_no);
            
            

            // for(int k=0;k<population_num;k++){
            //     if(cur_sol_fit[k]>=previous_sol_fit[k]){
                    
            //     }
            // }
            
            
            
            iter_obj_avg[j]+=best_score;
        }
        if(best_score<each_run_best_obj){
            each_run_best_obj=best_score;
            each_run_best=best_score_pos;
        }
        // cout<<"run:"<<i<<" "<<best_score<<endl;
    }
    ofstream output("result/mfo/"+func_name+".txt");
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