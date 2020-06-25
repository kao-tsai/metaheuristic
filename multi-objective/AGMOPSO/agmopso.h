#ifndef __AGMOPSO_H_INCLUDED__
#define __AGMOPSO_H_INCLUDED__
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

class agmopso{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
        typedef vector<int> d1;
        typedef vector<d1> d2;
		agmopso(int,int,int,int,double,double,string);
		void run();
	private:
		void init();
        int C(int,int);
		solution crossover(solution&,solution&);
        void mutation(solution&);
        population fitness(population&);
        population immune_search();
        void cal_pop_num();
        void generate_weight_vector();
        void get_near_T();
        void obj_sorting(solution&,vector<int>&);
        void update_z();
        double PBI(solution&,solution&);
        double VectorNorm2(solution&);
        double InnerProduct(solution&,solution&);
        population clonal_selection();
        void generate_z();
        void update_neighbor_sol(int,solution&,solution&);
        double TCH(solution&,solution&);
        double cal_euDis(solution&,solution&);
        int A_dominate_B(solution&,solution&);
        void archive_init();
        solution crowding_degree_cal();
        void quick_sort_obj(int,d1&,d1&,int,int);
        void quick_sort(solution&,d1&,int,int);
        void archive_update(population&);
        int select_pbest(int);
        int select_lbest(int);
        int select_gbest();
        solution PSO_search(int,int,int,int,solution&);
        double get_F1(int,int);
private:
	int numRuns;
	// int numIter;
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
    population archive;
    population archive_fit;
    int NC;
    int max_fes;
    int cur_fes;
    
};
agmopso::agmopso(int xNumRuns,
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
	max_fes=xNumIter;
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
int agmopso::C(int up,int down){
    if(down==0)
        return 1;
    else if(down==1)
        return up;
    else if(up==down)
        return 1;
    else
        return C(up-1,down-1)+C(up-1,down);
    
}
void agmopso::cal_pop_num(){
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
    if(obj_num==3)
        population_num=105;
    if(obj_num==2)
        population_num=100;
    // population_num=C(H+obj_num-1,obj_num-1);

}
void agmopso::init(){
    
    archive.clear();
    T_set.clear();
    ref_z.clear();
    nadir_z.clear();
    cur_fes=0;
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
    NC=population_num/5;
    archive_init();
}

void agmopso::archive_init(){
    archive=current_sol;
    archive_fit=current_sol_fit;
    //higer-objective
    d2 set(population_num,d1(1,0));
    
    for(int i=0;i<population_num-1;i++) 
        for(int j=i+1;j<population_num;j++){
            int do_count=0;
            int equal_count=0;
            for(int ob=0;ob<obj_num;ob++)
                if(archive_fit[i][ob]<=archive_fit[j][ob]){
                    do_count++;
                    if(archive_fit[i][ob]==archive_fit[j][ob])
                        equal_count++;
                }
            if(equal_count==obj_num)
                continue;
            else
                if(do_count==obj_num){
                    set[i].push_back(j);
                    set[j][0]++;
                }
                else if(do_count==0){
                    set[i][0]++;
                    set[j].push_back(i);
                }
        }
    //開始排rank
    d2 rank;
    d1 rank_i;
    bool not_empty;
    while(1){
        not_empty=false;
        for(int i=0;i<population_num;i++)       
            if(set[i][0]==0){
                rank_i.push_back(i);
                set[i][0]=-1;
                not_empty=true;
            }
        if(!not_empty)
            break;
        rank.push_back(rank_i);
        for(int i=0;i<rank_i.size();i++)
            for(int j=1;j<set[rank_i[i]].size();j++)
                set[set[rank_i[i]][j]][0]--;
        rank_i.clear();
    }
    population new_archive(rank[0].size());
    population new_archive_fit(rank[0].size());
    for(int i=0;i<rank[0].size();i++){
        new_archive[i]=archive[rank[0][i]];
        new_archive_fit[i]=archive_fit[rank[0][i]];
    }
    archive=new_archive;
    archive_fit=new_archive_fit;
}
//Return value 1 => A dominate B , 2 => B dominate A , 3 => A and B are equal , 4 => A and B are non-dominated 
int agmopso::A_dominate_B(solution& A_fit,solution& B_fit){
    int do_count=0;
    int equal_count=0;
    for(int ob=0;ob<obj_num;ob++)
        if(A_fit[ob]<=B_fit[ob]){
            do_count++;
            if(A_fit[ob]==B_fit[ob])
                equal_count++;
        }
    if(equal_count==obj_num)
        return 3;
    else
        if(do_count==obj_num)
            return 1;
        else if(do_count==0)
            return 2;
        else
            return 4;
        
}
void agmopso::generate_weight_vector(){
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

void agmopso::get_near_T(){
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
void agmopso::generate_z(){
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

void agmopso::update_z(){
    int archive_size=archive.size();
    for(int i=0;i<archive_size;i++)
        for(int j=0;j<obj_num;j++)
            if(archive_fit[i][j]<ref_z[j])
                ref_z[j]=archive_fit[i][j];
}

void agmopso::quick_sort(solution& value,vector<int>& this_rank,int left,int right){
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
double agmopso::cal_euDis(solution& t1,solution& t2){
    double sum=0;
    double len=t1.size();
    for(int i=0;i<len;i++)
        sum+=pow(t2[i]-t1[i],2.0);
    return sqrt(sum);
}
//-----------------Simulated Binary Crossover(From NSGAII SBX)------------------------//
agmopso::solution agmopso::crossover(solution& x1,solution& x2){
    
    random_device rd;
    mt19937 generator(rd());
    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1);
    // cout<<rndF(gen)<<endl;
    // cout<<rndF(gen)<<endl;
    // exit(0);
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

void agmopso::mutation(solution& x1){
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

agmopso::population agmopso::immune_search(){
    population parent;
    population child;
    
    parent = clonal_selection();
    child.assign(parent.size(),solution());
    int index1,index2;
    for(int i=0;i<parent.size();i++){
        index1=i;
        do{
            index2=rand()%population_num;
        }while(index1==index2);
        
        child[i]=crossover(parent[index1],parent[index2]);
        
        mutation(child[i]);
    }
    return child;
}
void agmopso::update_neighbor_sol(int cur_i,solution& new_sol,solution& new_sol_fit){
    double old_gtch,new_gtch;
    // cout<<cur_i<<endl;
    // cout<<weight_vector[T_set[cur_i][0]].size()<<endl;
    // exit(0);
    for(int i=0;i<T_num;i++){
        // new_gtch=TCH(new_sol_fit,weight_vector[T_set[cur_i][i]]);
        // old_gtch=TCH(current_sol_fit[T_set[cur_i][i]],weight_vector[T_set[cur_i][i]]);
    
        new_gtch=PBI(new_sol_fit,weight_vector[T_set[cur_i][i]]);
        old_gtch=PBI(current_sol_fit[T_set[cur_i][i]],weight_vector[T_set[cur_i][i]]);

        if(new_gtch<old_gtch){
            current_sol[T_set[cur_i][i]]=new_sol;
            current_sol_fit[T_set[cur_i][i]]=new_sol_fit;
        }
    }
}
double agmopso::TCH(solution& sol_fit,solution& cur_wei){
    //Without Normalized
    double max_num=fabs(sol_fit[0]-ref_z[0])*cur_wei[0];
    double obj_tch;
    for(int i=1;i<obj_num;i++){
        obj_tch=fabs(sol_fit[i]-ref_z[i])*cur_wei[i];
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
double agmopso::PBI(solution& sol_fit,solution& cur_wei){
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
    // double scalar_obj;
    // vector<double> vect_normal = vector<double>(obj_num,0);
	// 	double namda_norm = this->VectorNorm2(cur_wei);
	// 	for(int n=0; n<obj_num; n++)
	// 	{
	// 		vect_normal[n] = cur_wei[n]/namda_norm;
	// 	}			

	// 	double d1_proj = abs(this->InnerProduct(sol_fit, vect_normal) - this->InnerProduct(ref_z, vect_normal));
	// 	double temp = d1_proj*d1_proj-2*d1_proj*(this->InnerProduct(sol_fit,vect_normal) - this->InnerProduct(ref_z,vect_normal))
	// 		          + (this->InnerProduct(sol_fit, sol_fit) + this->InnerProduct(ref_z, ref_z) - 2*this->InnerProduct(sol_fit, ref_z));				
	// 	double d2_pred  = sqrt(temp);
	// 	scalar_obj = d1_proj + 20*d2_pred;
    // return scalar_obj;
}
double agmopso::InnerProduct(solution &vec1,solution &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=vec1[n]*vec2[n];
	return sum;
}
double agmopso::VectorNorm2(solution &vec1)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=vec1[n]*vec1[n];
	return sqrt(sum);
}
agmopso::population agmopso::clonal_selection(){
    population PC(NC);
    population PC_fit(NC);
    int archive_size=archive.size();
    solution crowding_value;
    d1 crowding_id(archive_size);

    for(int i=0;i<archive_size;i++)
        crowding_id[i]=i;


    // cout<<archive_size<<endl;
    if(archive_size>NC){
        PC.assign(NC,solution());
        PC_fit.assign(NC,solution());
        crowding_value=crowding_degree_cal();   
        quick_sort(crowding_value,crowding_id,0,archive_size-1);
        reverse(crowding_value.begin(),crowding_value.end());
        reverse(crowding_id.begin(),crowding_id.end());
        for(int i=0;i<NC;i++){
            PC[i]=archive[crowding_id[i]];
            PC_fit[i]=archive_fit[crowding_id[i]];
        }
    }
    else
    {
        
        PC.assign(archive_size,solution());
        PC_fit.assign(archive_size,solution());
        crowding_value=crowding_degree_cal();
        for(int i=0;i<archive_size;i++){
            PC[i]=archive[crowding_id[i]];
            PC_fit[i]=archive_fit[crowding_id[i]];
        }
    }
    
    int temp_size;
    if(archive_size>NC)
        temp_size=NC;
    else
        temp_size=archive_size;
    
    population new_pop_E;
    double crowding_sum=0.0;
    for(int i=0;i<temp_size;i++)
        crowding_sum+=crowding_value[i];
    
    // cout<<"this"<<crowding_value.size()<<endl;
    int q_i;
    
    int test_sum=0;
    for(int i=0;i<temp_size;i++){
        // cout<<crowding_value[i]/crowding_sum<<endl;
        q_i=ceil(((double)population_num)*(crowding_value[i]/crowding_sum));
        // cout<<i<<" "<<q_i<<endl;
        // test_sum+=q_i;
        for(int j=0;j<q_i;j++){
            new_pop_E.push_back(PC[i]);
            
        }
        
    }
    // cout<<"new sum:"<<test_sum<<endl;
    //-------------記得刪----------------//
    // for(int i=0;i<population_num;i++)
    // {
    //     for(int j=0;j<dimension;j++)
    //         cout<<new_pop_E[i][j]<<" ";
    //     cout<<endl;
    // }
    // exit(0);
    
    return new_pop_E;
}

agmopso::solution agmopso::crowding_degree_cal(){

    int front_size=archive_fit.size();
    solution init_dis(front_size,0.0);
    d1 this_rank(front_size);
    for(int i=0;i<front_size;i++)
        this_rank[i]=i;
    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    for(int i=0;i<obj_num;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    }
    if(front_size>2){
        //計算擁擠程度
        for (int i=0; i<obj_num; i++)
            init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
            
        for(int i=0;i<obj_num;i++)
            for(int j=1;j<front_size-1;j++){
                // if(init_dis[j]!=numeric_limits<double>::infinity()){
                    // if(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]== combine_fit[this_rank[sorting_id[i][0]]][i])
                    //     init_dis[sorting_id[i][j]]+=0.0;
                    // else
                        init_dis[sorting_id[i][j]]+=(archive_fit[this_rank[sorting_id[i][j+1]]][i]-archive_fit[this_rank[sorting_id[i][j-1]]][i])/
                        (archive_fit[this_rank[sorting_id[i][front_size-1]]][i]-archive_fit[this_rank[sorting_id[i][0]]][i]);
                // }
            }
        //距離密度值進行正規化
        for (int j=0; j<front_size; j++){
            if (init_dis[j]!= numeric_limits<double>::infinity())
            init_dis[j] = init_dis[j]/obj_num;
        }
        
        //two_times infinity
        double max_value=-numeric_limits<double>::infinity();
        for(int i=0;i<front_size;i++)
            if(init_dis[i]>max_value && init_dis[i]!=numeric_limits<double>::infinity())
                max_value=init_dis[i];
        for(int i=0;i<front_size;i++)
            if(init_dis[i]==numeric_limits<double>::infinity())
                init_dis[i]=max_value*2;

        // for (int j=0; j<front_size; j++){
        //     cout<<init_dis[j]<<endl;;
        // }
        // exit(0);
    }
    //when front size = 2
    else if(front_size==2){
        init_dis[0]=0.5;
        init_dis[1]=0.5;
    }
    else if(front_size==1)
        init_dis[0]=1;

    return init_dis;
}

void agmopso::quick_sort_obj(int obj_level,d1& obj_id,d1& sorting_id,int left,int right){
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        temp = obj_id[right];
        obj_id[right] = obj_id[index];
        obj_id[index] = temp;

        temp = sorting_id[right];
        sorting_id[right] = sorting_id[index];
        sorting_id[index] = temp;

        pivot = archive_fit[obj_id[right]][obj_level];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (archive_fit[obj_id[j]][obj_level] <= pivot)
            {
                i+=1;
                temp = obj_id[j];
                obj_id[j] = obj_id[i];
                obj_id[i] = temp;
                
                temp = sorting_id[j];
                sorting_id[j] = sorting_id[i];
                sorting_id[i] = temp;
            }
        }
        index=i+1;
        temp = obj_id[index];
        obj_id[index] = obj_id[right];
        obj_id[right] = temp;

        temp = sorting_id[index];
        sorting_id[index] = sorting_id[right];
        sorting_id[right] = temp;

        quick_sort_obj(obj_level,obj_id,sorting_id, left, index-1);
        quick_sort_obj(obj_level,obj_id,sorting_id, index+1, right);
    }
}
void agmopso::archive_update(population& new_sol){
    
    int new_sol_size=new_sol.size();
    // if(cur_fes>0)
    //     cout<<new_sol_size<<endl;
    population new_sol_fit=fitness(new_sol);
    // if(cur_fes>0)
    //     cout<<"hello"<<endl;
    int archive_size;
    population temp_archive,temp_archive_fit;
    d1 mark(archive_size,0);
    int flag;
    
    for(int i=0;i<new_sol_size;i++){
        archive_size=archive.size();
        mark.assign(archive_size,0);
        temp_archive=archive;
        temp_archive_fit=archive_fit;
        for(int j=0;j<archive_size;j++){
            flag=A_dominate_B(new_sol_fit[i],archive_fit[j]);
            if(flag==1)
                mark[j]=1;
            else if(flag==2)
                break;
        }
        
        archive.clear();
        archive_fit.clear();
        for(int k=0;k<mark.size();k++)
            if(mark[k]==0){
                archive.push_back(temp_archive[k]);
                archive_fit.push_back(temp_archive_fit[k]);
            }
        solution cur_crowding_value;
        if(flag!=2){
            archive.push_back(new_sol[i]);
            archive_fit.push_back(new_sol_fit[i]);
            
            if(archive.size()>population_num){
                cur_crowding_value=crowding_degree_cal();
                double most_crowded_value=cur_crowding_value[0];
                int most_crowded_id=0;
                for(int i=0;i<archive.size();i++)
                    if(most_crowded_value>cur_crowding_value[i]){
                        most_crowded_value=cur_crowding_value[i];
                        most_crowded_id=i;
                    }
                archive.erase(archive.begin()+most_crowded_id);
                archive_fit.erase(archive_fit.begin()+most_crowded_id);
            }
        }
    }

}

int agmopso::select_pbest(int cur_j){
    int min_pbest_j;

    //Select personal best
    min_pbest_j=cur_j;
    for(int j=0;j<archive.size();j++)
        if(PBI(archive_fit[min_pbest_j],weight_vector[cur_j])>PBI(archive_fit[j],weight_vector[cur_j]))
            min_pbest_j=j;
    return min_pbest_j;
}
int agmopso::select_lbest(int cur_j){
    int near_id;
    near_id=rand()%T_num;
    int min_lbest_j;
    min_lbest_j=cur_j;
    for(int j=0;j<archive.size();j++)
        if(PBI(archive_fit[min_lbest_j],weight_vector[T_set[cur_j][near_id]])>PBI(archive_fit[j],weight_vector[T_set[cur_j][near_id]]))
            min_lbest_j=j;
    return min_lbest_j;
}
int agmopso::select_gbest(){
    int r=rand()%archive.size();
    return r;
}
agmopso::solution agmopso::PSO_search(int cur_j,int pbest_id,int lbest_id,int gbest_id,solution& Vij){
    double F1=get_F1(cur_j,pbest_id);
    double F2=0.5;
    double w=0.1;
    w=((double)rand()/RAND_MAX)*0.4+0.1;
    solution new_sol(dimension);
    for(int i=0;i<dimension;i++){
        Vij[i]=w*Vij[i]+F1*(archive[pbest_id][i]-current_sol[cur_j][i])+F2*(archive[lbest_id][i]-archive[gbest_id][i]);
        new_sol[i]=current_sol[cur_j][i]+Vij[i];
        //Keep solution in the bound
        if(new_sol[i]>upperbound[i])
            new_sol[i]=upperbound[i];
        else if(new_sol[i]<lowerbound[i])
            new_sol[i]=lowerbound[i];
    }
    return new_sol;
}
double agmopso::get_F1(int cur_j,int pbest_id){
    double dis1;
    double tmp_sum=0,delta=20;
    double norm_vec;
    solution tmp_zero(obj_num,0.0);
    for(int i=0;i<obj_num;i++)
        tmp_sum+=fabs((archive_fit[pbest_id][i]-ref_z[i])*weight_vector[cur_j][i]);
    norm_vec=cal_euDis(weight_vector[cur_j],tmp_zero);
    dis1=tmp_sum/norm_vec;
    return dis1;
}
//-----------------------------------------------------run MOEA/D algorithm----------------------------------------------------------//
void agmopso::run(){
    
    population S;
    population P;
    
    population temp_2d(1);
    population temp_fit;
    
    int pbest_id;
    int lbest_id;
    int gbest_id;
    population vel;
    for(int i=0;i<numRuns;i++){
        
        init();
        P.assign(population_num,solution());
        vel=population(population_num,solution(dimension,0));
        while(cur_fes<max_fes){
            
            S=immune_search();
            // cout<<S.size()<<endl;
            
            archive_update(S);
            
            cur_fes+=S.size();
            update_z();
            
            for(int j=0;j<population_num;j++){
                
                pbest_id=select_pbest(j);
                lbest_id=select_lbest(j);
                gbest_id=select_gbest();
                // cout<<j<<" hello"<<endl;
                P[j]=PSO_search(j,pbest_id,lbest_id,gbest_id,vel[j]);
                
            }
            // cout<<P.size()<<endl;
            
            archive_update(P);
            
            cur_fes+=P.size();
            current_sol=P;
            update_z();

        }
        
       
        // //輸出該RUN的所有解
        ofstream output_obj;
        output_obj.open("pareto/"+func_name+"/agmopso/"+func_name+"_agmopso_"+to_string(i)+".txt");
          for(int k=0;k<archive_fit.size();k++){
            for(int l=0;l<obj_num;l++)
                output_obj<<archive_fit[k][l]<<" ";
            output_obj<<endl;
        }
        output_obj.close();
    }
    
   

}

agmopso::population agmopso::fitness(population& all_sol) {
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
    else if(func_name=="UF2")
        return func.UF2(all_sol);
    else if(func_name=="UF3")
        return func.UF3(all_sol);
    else if(func_name=="UF4")
        return func.UF4(all_sol);
    else if(func_name=="UF5")
        return func.UF5(all_sol);
    else if(func_name=="UF6")
        return func.UF6(all_sol);
    else if(func_name=="UF7")
        return func.UF7(all_sol);
    else if(func_name=="UF8")
        return func.UF8(all_sol);
    else if(func_name=="UF9")
        return func.UF9(all_sol);
    else if(func_name=="UF10")
        return func.UF10(all_sol);
    else if(func_name=="DTLZ1")
        return func.DTLZ1(all_sol);
    else if(func_name=="DTLZ2")
        return func.DTLZ2(all_sol);
    else if(func_name=="DTLZ3")
        return func.DTLZ3(all_sol);
    else if(func_name=="DTLZ4")
        return func.DTLZ4(all_sol);
    else if(func_name=="DTLZ5")
        return func.DTLZ5(all_sol);
    else if(func_name=="DTLZ6")
        return func.DTLZ6(all_sol);
    else if(func_name=="DTLZ7")
        return func.DTLZ7(all_sol);
}




#endif
       