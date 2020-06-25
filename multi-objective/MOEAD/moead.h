#ifndef __MOEAD_H_INCLUDED__
#define __MOEAD_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <random>
#include "../test_problem/test_problem.h"

using namespace std;
# define EPS 1.0e-14

class moead{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		moead(int,int,int,int,double,double,string);
		void run();
	private:
		void init(int);
        int C(int,int);
		solution crossover(solution&,solution&);
        void mutation(solution&);
        population fitness(population&);
        solution transition(int);
        void cal_pop_num();
        void generate_weight_vector();
        void get_near_T();
        void obj_sorting(solution&,vector<int>&);
        void quick_sort(vector<double>&,vector<int>&,int,int);
        void update_z(int,solution&);
        double PBI(solution&,solution&);
        // double VectorNorm2(solution&);
        // double InnerProduct(solution&,solution&);
        void generate_z();

        void update_neighbor_sol(int,solution&,solution&);
        double TCH(solution&,solution&);
        double cal_euDis(solution&,solution&);
private:
	int numRuns;
	int numIter;
    int population_num;
	double mutation_pro;
    double crossover_pro;
    string func_name;
    int dimension; //每一個解的維度
    int obj_num; //目標數量
    solution upperbound;
    solution lowerbound;
    double dimm;
    double dimc;
    vector<int> nbits;
    solution ref_z;
    solution nadir_z;
    population weight_vector;
    int H;
    int T_num;
    vector<vector<int>> T_set;
    test_problem func;

    population current_sol;
    population current_sol_fit;
    double inf=numeric_limits<double>::infinity();
    
    
};
moead::moead(int xNumRuns,
int xNumIter,
int xH,
int xT_num,
double xcrossover_pro,
double xmutation_pro,
string xfunc_name
)
{
    
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
    crossover_pro=xcrossover_pro;
    // mutation_pro=xmutation_pro;
    func_name=xfunc_name;
    H=xH;
    T_num=xT_num;
    func=test_problem(func_name);
    dimension=func.idimension;
    obj_num=func.iobj_num;
    upperbound=func.ub;
    lowerbound=func.lb;
    mutation_pro=1.0/(double)dimension;
    
    dimm=20;//real-coded ngsaii
    dimc=20;//real-coded nsgaii
    
}
// int moead::C(int up,int down){
//     if(down==0)
//         return 1;
//     else if(down==1)
//         return up;
//     else if(up==down)
//         return 1;
//     else
//         return C(up-1,down-1)+C(up-1,down);
    
// }
void moead::cal_pop_num(){
    // int HM1=1,HM2=1,HM3=1;
    // for(int i=1;i<=H+obj_num-1;i++){
    //     if(i<=H+obj_num-1){
    //         HM1*=i;
    //         cout<<HM1<<endl;
    //     }
    //     if(i<=obj_num-1)
    //         HM2*=i;
    //     if(i<=H)
    //         HM3*=i;
    // }
    
    population_num=300;
    // population_num=C(H+obj_num-1,obj_num-1);

}

void moead::init(int cur_run){
    current_sol.clear();
    current_sol_fit.clear();
    weight_vector.clear();
    T_set.clear();
    ref_z.clear();
    nadir_z.clear();
    // if(cur_run==0){
        cal_pop_num();//Calculate the population number with H
        
        current_sol=population(population_num,solution(dimension,0.0));
        current_sol_fit=population(population_num,solution(obj_num,0.0));
        weight_vector.assign(population_num,solution(obj_num));
    
        generate_weight_vector();
        
        T_set.assign(population_num,vector<int>(T_num));
        get_near_T();
    
    // }
	for (int i = 0; i < population_num; i++)
        for(int j=0;j<dimension;j++)
            current_sol[i][j]=(upperbound[j]-lowerbound[j])*((double)rand()/RAND_MAX)+lowerbound[j];

    current_sol_fit=fitness(current_sol);
    
    ref_z.assign(obj_num,0);
    nadir_z.assign(obj_num,0);
    generate_z();//Generate the ideal point

    
}
void moead::generate_weight_vector(){
    solution plus_pos(obj_num-1,0);
    solution current_set(obj_num);
    int final_pos=plus_pos.size()-1;
    int pos1=0,pos2=H;
    int k=0,accept_pos;
    while(plus_pos[0]<=H){
        for(int i=0;i<=final_pos;i++){
            if(i==0)
                current_set[i]=(plus_pos[i]-pos1)/H;
            else
                current_set[i]=(plus_pos[i]-plus_pos[i-1])/H;
            
            if(i==final_pos)
                current_set[i+1]=(pos2-plus_pos[i])/H;
        }
        weight_vector[k]=current_set;
        k++;
        //get new set
        plus_pos[final_pos]++;
        for(int i=0;i<obj_num-1;i++)
            if(plus_pos[final_pos-i]<=H){
                accept_pos=final_pos-i;
                break;
            }
            else if(plus_pos[final_pos-i]>H && final_pos-i!=0)
                plus_pos[final_pos-i-1]++;
            else{
                accept_pos=final_pos;
                break;
            }
        
        for(int i=accept_pos+1;i<=final_pos;i++)
            plus_pos[i]=plus_pos[accept_pos];
        
    }
}

void moead::get_near_T(){
    vector<vector<int>> sorting_id(population_num,vector<int>(population_num));
    for(int i=0;i<population_num;i++)
        for(int j=0;j<population_num;j++)
            sorting_id[i][j]=j;
    
    population weight_vector_dis(population_num,solution(population_num));
    double tmp_dis;
    for(int i=0;i<population_num;i++)
        for(int j=i;j<population_num;j++){
            tmp_dis=cal_euDis(weight_vector[i],weight_vector[j]);
            weight_vector_dis[i][j]=tmp_dis;
            weight_vector_dis[j][i]=tmp_dis;
        }
    
    for(int i=0;i<population_num;i++){
        quick_sort(weight_vector_dis[i],sorting_id[i],0,population_num-1);
        for(int j=0;j<T_num;j++){
            T_set[i][j]=sorting_id[i][j];
           
        }
        
    }
    
}
void moead::generate_z(){
    for(int i=0;i<obj_num;i++){
        ref_z[i]=inf;
        nadir_z[i]=-inf;
    }
    for(int i=0;i<population_num;i++)
        for(int j=0;j<obj_num;j++){
            if(current_sol_fit[i][j]<ref_z[j])
                ref_z[j]=current_sol_fit[i][j];
            if(current_sol_fit[i][j]>nadir_z[j])
                nadir_z[j]=current_sol_fit[i][j];
        }
   
}

void moead::update_z(int cur_i,solution& new_sol_fit){
    for(int i=0;i<obj_num;i++){
        nadir_z[i]=-inf;
        if(new_sol_fit[i]<ref_z[i])
            ref_z[i]=new_sol_fit[i];
    }
    for(int i=0;i<population_num;i++)
        for(int j=0;j<obj_num;j++){
            if(cur_i==i)
                if(new_sol_fit[j]>nadir_z[j])
                    nadir_z[j]=new_sol_fit[j];        
            if(current_sol_fit[i][j]>nadir_z[j])
                nadir_z[j]=current_sol_fit[i][j];
        }
}

void moead::quick_sort(solution& value,vector<int>& this_rank,int left,int right){
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
double moead::cal_euDis(solution& t1,solution& t2){
    double sum=0;
    double len=t1.size();
    for(int i=0;i<len;i++)
        sum+=pow(t2[i]-t1[i],2.0);
    return sqrt(sum);
}
//-----------------Simulated Binary Crossover(From NSGAII SBX)------------------------//
moead::solution moead::crossover(solution& x1,solution& x2){
    
    random_device rd;
    mt19937 generator(rd());
    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1);
    
    double rnd,betaq,beta,alpha;
    solution child(dimension);
    double c1,c2;
    double y1,y2,yu,yl;
    
    rnd =rndF(gen);
    /*Check Whether the cross-over to be performed*/
    if(rnd <= crossover_pro){
        
        /*Loop over no of variables*/
        for(int j = 0;j < dimension;j++){
        /*Selected Two Parents*/  
            rnd = rndF(gen);             
            /* Check whether variable is selected or not*/
            if(rnd <= 0.5){
                /*Variable selected*/
                if(fabs(x1[j] - x2[j]) > EPS){
                    if(x2[j] > x1[j]){
                        y2 = x2[j];
                        y1 = x1[j];
                    }
                    else{
                        y2 = x1[j];
                        y1 = x2[j];
                    }
                    yl = lowerbound[j];
                    yu = upperbound[j];
                    
                    rnd=rndF(gen);
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(dimc+1.0));
                    
                    if (rnd <= 1.0/alpha)
                    {
                        betaq = pow ((rnd*alpha),(1.0/(dimc+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rnd*alpha)),(1.0/(dimc+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(dimc+1.0));
                        if (rnd <= (1.0/alpha))
                    {
                        betaq = pow ((rnd*alpha),(1.0/(dimc+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rnd*alpha)),(1.0/(dimc+1.0)));
                    }
                        c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    rnd=rndF(gen);
                    if (rnd<=0.5)
                    {
                        child[j]=c1;
                    }
                    else
                    {
                        child[j]=c2;
                    }
                }
                else{
                    /*Copying the children to parents*/
                    child[j]=x1[j];
                }
            }
            else{
                child[j]=x1[j];
            }
        }

    }
    else
        for (int j=0; j<dimension; j++)     
            child[j] = x1[j];
        
    return child;

}

void moead::mutation(solution& x1){
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    random_device rd;
    mt19937 generator(rd());
    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1);
    for (j=0; j<dimension; j++)
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
        }
    }

}

moead::solution moead::transition(int cur_i){
    int index1,index2;
    index1=rand()%T_num;
    do{
        index2=rand()%T_num;
    }while(index1==index2);
    
    solution child;
    child=crossover(current_sol[T_set[cur_i][index1]],current_sol[T_set[cur_i][index2]]);
    mutation(child);
    
    return child;
}
void moead::update_neighbor_sol(int cur_i,solution& new_sol,solution& new_sol_fit){
    double old_gtch,new_gtch;
    
    for(int i=0;i<T_num;i++){

        new_gtch=TCH(new_sol_fit,weight_vector[T_set[cur_i][i]]);
        old_gtch=TCH(current_sol_fit[T_set[cur_i][i]],weight_vector[T_set[cur_i][i]]);
    
        // new_gtch=PBI(new_sol_fit,weight_vector[T_set[cur_i][i]]);
        // old_gtch=PBI(current_sol_fit[T_set[cur_i][i]],weight_vector[T_set[cur_i][i]]);

        if(new_gtch<old_gtch){
            current_sol[T_set[cur_i][i]]=new_sol;
            current_sol_fit[T_set[cur_i][i]]=new_sol_fit;
        }
    }
}
double moead::TCH(solution& sol_fit,solution& cur_wei){
    //Without Normalized
    // double max_num=fabs(sol_fit[0]-ref_z[0])*cur_wei[0];
    
    // double obj_tch;
    // for(int i=0;i<obj_num;i++){
    //     obj_tch=fabs(sol_fit[i]-ref_z[i])*cur_wei[i];
    //     if(obj_tch>max_num)
    //         max_num=obj_tch;
    // }
    // return max_num;
    //Modified TCH()
    double max_num=fabs(sol_fit[0]-ref_z[0])/cur_wei[0];
    
    double obj_tch;
    for(int i=0;i<obj_num;i++){
        obj_tch=fabs(sol_fit[i]-ref_z[i])/cur_wei[i];
        if(obj_tch>max_num)
            max_num=obj_tch;
    }
    return max_num;
    //Normalized
    // double max_num=fabs((sol_fit[0]-ref_z[0])/(nadir_z[0]-ref_z[0]))*cur_wei[0];
    // double obj_tch;
    // for(int i=1;i<obj_num;i++){
    //     obj_tch=fabs((sol_fit[i]-ref_z[i])/(nadir_z[i]-ref_z[i]))*cur_wei[i];
    //     if(obj_tch>max_num)
    //         max_num=obj_tch;
    // }
    // return max_num;
}
double moead::PBI(solution& sol_fit,solution& cur_wei){
    double dis1,dis2;
    double tmp_sum=0,delta=20;
    double norm_vec;
    solution tmp_zero(obj_num,0.0);
    for(int i=0;i<obj_num;i++)
        tmp_sum+=fabs((sol_fit[i]-ref_z[i])*cur_wei[i]);
    norm_vec=cal_euDis(cur_wei,tmp_zero);
    dis1=tmp_sum/norm_vec;
    tmp_sum=0;
    solution dis2_vec(obj_num,0.0);
    for(int i=0;i<obj_num;i++)
        dis2_vec[i]=(sol_fit[i]-(ref_z[i]+dis1*(cur_wei[i]/norm_vec)));
    dis2=cal_euDis(dis2_vec,tmp_zero);
    return dis1+dis2*delta;
    
}


//-----------------------------------------------------run MOEA/D algorithm----------------------------------------------------------//
void moead::run(){
    
    solution child(dimension,0.0);
    solution child_fit(obj_num,0);
    population temp_2d(1);
    population temp_fit;
    
    for(int i=0;i<numRuns;i++){
        
        init(i);
        population output_archive(population_num);
        for(int j=0;j<numIter;j++){
            
            for(int k=0;k<population_num;k++){
                child=transition(k);
                
                // //-----------------------//
                temp_2d[0]=child;
                temp_fit=fitness(temp_2d);
                child_fit=temp_fit[0];
                
                // //-----------------------//
                update_z(k,child_fit);
                
                update_neighbor_sol(k,child,child_fit);
                
                
            }

        }
        
       
        // //輸出該RUN的所有解
        ofstream output_obj;
        output_obj.open("pareto/"+func_name+"/moead/"+func_name+"_moead_"+to_string(i)+".txt");
          for(int k=0;k<population_num;k++){
            for(int l=0;l<obj_num;l++)
                output_obj<<current_sol_fit[k][l]<<" ";
            output_obj<<endl;
        }
        output_obj.close();
    }
    
   

}

moead::population moead::fitness(population& all_sol) {
    if(func_name=="SCH")
        return func.SCH(all_sol);
    else if(func_name=="FON")
        return func.FON(all_sol);
    else if(func_name=="POL")
        return func.POL(all_sol);
    else if(func_name=="KUR")
        return func.KUR(all_sol);
    else if(func_name=="ZDT1")
        return func.ZDT1(all_sol);
    else if(func_name=="ZDT2")
        return func.ZDT2(all_sol);
    else if(func_name=="ZDT3")
        return func.ZDT3(all_sol);
    else if(func_name=="ZDT4")
        return func.ZDT4(all_sol);
    else if(func_name=="ZDT6")
        return func.ZDT6(all_sol);
    else if(func_name=="UF1")
        return func.UF1(all_sol);
    else if(func_name=="DTLZ1")
        return func.DTLZ1(all_sol);
    else if(func_name=="DTLZ2")
        return func.DTLZ2(all_sol);
}




#endif
       