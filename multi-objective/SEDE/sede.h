#ifndef __SEDE_H_INCLUDED__
#define __SEDE_H_INCLUDED__
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

class sede{
	public:
		typedef vector<int> d1;
        typedef vector<d1> d2;
        typedef vector<d2> d3;
        typedef vector<d3> d4;
        typedef vector<double> dd1;
        typedef vector<dd1> dd2;
        typedef vector<dd2> dd3;
        typedef vector<dd3> dd4;
		sede(int,int,int,int,int,int,int,string,double,double);
		void run();
	private:
		void init();
        void ga_crossover();
        void de_mut();
        void spiral_tran();
        void create_rotMat();
        void arrange();
        void cal_ev();
        void choose_sample(int);
        double bound_rand(double,double);
        void select_player();
        void marketing_survey();
        void divide_region();
        dd1 decode(d2&);
        void cal_cdr_fit();
        dd2 fitness(dd2&);
        void crowding_dis_assign(d1&);
        void crowding_entropy(d1&);
        void Euclidean_crowding(d1&);
        double Euclidean_dis(dd1&,dd1&);
        void quick_sort(dd1&,d1&,int,int);
        void quick_sort_obj(int,d1&,d1&,int,int);
        void archive_delete();
        void crowding_d_evaluate(d1&);
        void Euclidean_archive_d(d1&);
        void Euclidean_d_evaluate(d1&);
        void combine_region();
        void update_uCr_uF();
        void generate_new_Cr_F();
        //Decomposiotion based
        void decom_archive_delete();
        void generate_weight_vector();
        void generate_z(dd2&);
        void update_z(dd2&);
        double TCH(dd1&,dd1&);
        dd1 return_crowd_val(d1&);
        void archive_update(dd2&);
        //
        int return_weight_belong(dd1&);
private:
	int numRuns;
	int numIter;
	int dim;

    dd2 searcher;
    dd3 sample;
    dd4 sampleV;
    d2 region_bit;//切法(1).針對維度少的function,
    dd3 region_bound;//切法(2).針對維度多的function
    dd1 ta;
    dd1 tb;
    d1 searcher_belong;
    dd2 EV;
    dd1 region_best_fit;//考慮是否改成cdr(當前最佳非支配解)
    
    // dd2 sample_fit;//各目標的Objective value()
    // dd2 sampleV_fit;//各目標的Objective value()
    // dd2 searcher_fit;//各目標的Objective value()
    
    d1 nbits;
    
    dd1 upperbound;
    dd1 lowerbound;

    dd2 archive;//Best Non-Dominated Solutions
    d1 archive_region_id;
    dd2 combine_sol;//暫存Sample+SampleV+Archive
    dd2 combine_fit;
    dd1 combine_cdr_fit;
    d2 combine_region_id_2d;
    d1 combine_region_id;
    d1 rank_0;
    int sampleV_start;
    dd2 archive_fit;//各目標的Objective value
    // dd2 archive_cdr_fit;//sample+sampleV+archive(Region,CDR)

    dd2 sample_cdr_fit;//參考Crowding Distance和Rank合併的Fitness(用來計算期望值公式)
    dd3 sampleV_cdr_fit;//參考Crowding Distance和Rank合併的Fitness(用來計算期望值公式)(Searcher,Region,CDR)
    d1  region_selected_by_s;
    int origin_region_num;
    int origin_sample_num;
    int origin_searcher_num;
    int full_sample_num;
    int max_pareto;
    int obj_num;
    int best_obj_val;
    int clip_bit_num;
    int region_num;
    int searcher_num;
    int sample_num;
    int tour_size;
    double crossover_pro;
    double mutation_pro;
    double tatb_pressure;
    double dominate_pressure;
    double crowd_p;
    double crowd_p_max;
    dd2 CR_i;
    dd2 F_i;
    dd1 S_cr;
    dd1 S_f;
    double u_cr;
    double u_f;
    int tmp_iter;
    string func_name;
    test_problem func;
    int deleted_pos;
    double angle;
    double r_dis;
    dd2 rotMat;
    dd2 rotMat_I;
    double inf=numeric_limits<double>::infinity();
    //Decomposition Based Archive Update
    dd2 weight_vector;
    dd1 ref_z;
    dd1 nadir_z;
    int H;
    int weight_vector_size;
    // d1 decom_region_id;
};

sede::sede(int xNumRuns,
int xNumIter,
int in_max_pareto,
int xRegion,
int in_searcher_num,
int in_sample_num,
int in_tournament_size,
string in_func_name,
double in_crossover_pro,
double in_mutation_pro
)
{
    srand(time(0));

    numRuns=xNumRuns;
	numIter=xNumIter;
    origin_region_num=xRegion;
    region_num=xRegion;
    searcher_num=in_searcher_num;
    origin_searcher_num=in_searcher_num;
    sample_num=in_sample_num;
    origin_sample_num=in_sample_num;
    tour_size=in_tournament_size;
    func_name=in_func_name;
    func=test_problem(func_name);
    dim=func.idimension;
    obj_num=func.iobj_num;
    max_pareto=in_max_pareto;
    lowerbound=func.lb;
    upperbound=func.ub;
    full_sample_num=sample_num*region_num;
    // tatb_pressure=2;
    dominate_pressure=3;//3
    crossover_pro=in_crossover_pro;
    mutation_pro=in_mutation_pro;
    // mutation_pro=1.0/(nbits[0]*dim);
    angle=M_PI/2.0;
    r_dis=0.85;
    if(obj_num==2)
        H=99;
    else if(obj_num==3)
        H=13;

}

void sede::generate_weight_vector(){
    dd1 plus_pos(obj_num-1,0);
    dd1 current_set(obj_num);
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
    weight_vector_size=weight_vector.size();
}
void sede::generate_z(dd2& all_sol_fit){
    for(int i=0;i<obj_num;i++){
        ref_z[i]=inf;
        nadir_z[i]=-inf;
    }
    //search searcher ideal point
    for(int i=0;i<all_sol_fit.size();i++)
        for(int j=0;j<obj_num;j++){
            if(all_sol_fit[i][j]<ref_z[j])
                ref_z[j]=all_sol_fit[i][j];
            if(all_sol_fit[i][j]>nadir_z[j])
                nadir_z[j]=all_sol_fit[i][j];
        }
}
void sede::update_z(dd2& new_sol_fit){
    for(int i=0;i<new_sol_fit.size();i++){
        for(int j=0;j<obj_num;j++){
            nadir_z[j]=-inf;
            if(new_sol_fit[i][j]<ref_z[j])
                ref_z[j]=new_sol_fit[i][j];
        }
    }
    for(int i=0;i<new_sol_fit.size();i++)
        for(int j=0;j<obj_num;j++){
            if(new_sol_fit[i][j]>nadir_z[j])
                nadir_z[j]=new_sol_fit[i][j];        
        }
}
//修改
void sede::init(){
    region_num=origin_region_num;
    sample_num=origin_sample_num;
    searcher_num=origin_searcher_num;
    tb.assign(region_num,0.0);
    ta.assign(region_num,1.0);
    searcher_belong.assign(searcher_num,0);
    searcher.assign(searcher_num,dd1(dim,0.0));
    sample.assign(region_num,dd2(sample_num,dd1(dim,0.0)));
    sampleV.assign(searcher_num,dd3(region_num,dd2(sample_num,dd1(dim,0.0))));
    sample_cdr_fit.assign(region_num,dd1(sample_num,0.0));
    sampleV_cdr_fit.assign(searcher_num,dd2(region_num,dd1(sample_num,0.0)));
    region_best_fit.assign(region_num,0);
    combine_region_id.assign(searcher_num+sample_num*region_num,-1);
    archive_region_id.assign(searcher_num+sample_num*region_num,-1);
    region_selected_by_s.assign(searcher_num,-1);
    CR_i.assign(region_num,dd1(sample_num,0.2));//0.2
    F_i.assign(region_num,dd1(sample_num,0.8));//0.6
    S_cr.clear();
    S_f.clear();
    u_cr=0.2;
    u_f=0.9;
    weight_vector.assign(100,dd1(obj_num));//大小要設定
    ref_z.assign(obj_num,0.0);
    nadir_z.assign(obj_num,0.0);
    
    // generate_weight_vector();
    
    //擴充upperbound和lowerbound
    // for(int i=0;i<dim;i++){
    //     upperbound[i]+=upperbound[i]*0.5;
    //     lowerbound[i]+=lowerbound[i]*0.5;
    // }
    archive=dd2(searcher_num+region_num*sample_num);
    
	for (int i = 0; i < searcher_num; i++){
        for(int j=0;j<dim;j++){
            searcher[i][j]=bound_rand(lowerbound[j],upperbound[j]);
            // cout<<searcher[i][j]<<" ";        
        }
        // cout<<endl;
        archive[i]=searcher[i];
    }
    
    // create_rotMat();
    
}

double sede::bound_rand(double lw,double ub){
    return ((double)rand()/RAND_MAX)*(ub-lw)+lw;
}

//This divided method is for high dimension problem
void sede::divide_region(){
    //Lower Bound在上，Upper Bound在下
    region_bound.assign(region_num,dd2(2,dd1(clip_bit_num,0)));
    double divide_two;
    for(int i=0;i<region_num;i++){
        int bin=i;     
        for(int j=clip_bit_num-1;j>=0;j--){
            
            if(bin%2==0){
                divide_two=(upperbound[j]+lowerbound[j])/2;
                region_bound[i][0][j]=lowerbound[j];
                region_bound[i][1][j]=divide_two+(upperbound[j]-lowerbound[j])*0.1;//10
            }
            else if(bin%2==1){
                divide_two=(upperbound[j]+lowerbound[j])/2;
                region_bound[i][0][j]=divide_two-(upperbound[j]-lowerbound[j])*0.1;//10
                region_bound[i][1][j]=upperbound[j];
            }
            bin/=2;
        }
    }
}
//修改
void sede::arrange(){
    clip_bit_num=log2(region_num);
    
    divide_region();
    
    //--------------------分配searcher到各個Region---------------------//
    for(int i=0;i<searcher_num;i++){
        searcher_belong[i]=i%region_num;
        for(int j=0;j<clip_bit_num;j++)
            searcher[i][j]=bound_rand(region_bound[searcher_belong[i]][0][j],region_bound[searcher_belong[i]][1][j]);
        combine_region_id[i]=searcher_belong[i];
        region_selected_by_s[i]=searcher_belong[i];
    }
    
    //----------Generate Samples and Assign them to each Region----------//
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            for(int k=0;k<dim;k++)
                if(k<clip_bit_num){
                    sample[i][j][k]=bound_rand(region_bound[i][0][k],region_bound[i][1][k]);
                    // cout<<sample[i][j][k]<<" ";
                }
                else{
                    sample[i][j][k]=bound_rand(lowerbound[k],upperbound[k]);
                    // cout<<sample[i][j][k]<<" ";
                }
            // cout<<endl;
            archive[searcher_num+i*sample_num+j]=sample[i][j];
            combine_region_id[searcher_num+i*sample_num+j]=i;
        }
    for(int i=0;i<searcher_num;i++){
        tb[searcher_belong[i]]++;
        ta[searcher_belong[i]]=1.0;
    }
    
    EV.assign(searcher_num,dd1(region_num,0.0));
    
    //Calculate Objective Value of each solution
    int archive_size=archive.size();
    archive_fit.clear();//可能不需要
    archive_fit=fitness(archive);
    // generate_z(archive_fit);
    //紀錄dominated數量和dominate set
    //two-objective
    // d2 set(archive_size,d1(1,0));
    // for(int i=0;i<archive_size-1;i++) 
    //     for(int j=i+1;j<archive_size;j++)
    //         if(archive_fit[i][0]<=archive_fit[j][0]&&archive_fit[i][1]<=archive_fit[j][1])
    //             if(archive_fit[i][0]==archive_fit[j][0]&&archive_fit[i][1]==archive_fit[j][1])
    //                 continue;
    //             else{
    //                 set[i].push_back(j);
    //                 set[j][0]++;
    //             }
    //         else if(archive_fit[i][0]>=archive_fit[j][0]&&archive_fit[i][1]>=archive_fit[j][1]){
    //             set[i][0]++;
    //             set[j].push_back(i);
    //         }

    //higer-objective
    d2 set(archive_size,d1(1,0));
    
    for(int i=0;i<archive_size-1;i++) 
        for(int j=i+1;j<archive_size;j++){
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
        for(int i=0;i<archive_size;i++)       
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
    
    combine_sol=archive;
    combine_fit=fitness(combine_sol);
    //更新最小以及最大參考點
    update_z(combine_fit);
    //做crowding_distance
    for(int i=0;i<rank.size();i++)
        if(rank[i].size()>2){
            crowding_dis_assign(rank[i]);
            // crowding_entropy(rank[i]);
        }
    
    //計算所有解的CDR
    double giving_level=(double)archive_size;
    double rank_size=(double)rank.size();
    combine_cdr_fit.assign(archive_size,0.0);
    for(int i=0;i<rank_size;i++)
        for(int j=0;j<rank[i].size();j++){
            combine_cdr_fit[rank[i][j]]=giving_level*pow(dominate_pressure,(rank_size-i));
            giving_level=giving_level-1;
        }
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++)
            sample_cdr_fit[i][j]=combine_cdr_fit[searcher_num+i*sample_num+j];
    
    rank_0=rank[0];
   
    dd2 tmp_archive;
    d1 tmp_region_id;
    tmp_archive=archive;
    tmp_region_id=combine_region_id;
    archive=dd2(rank_0.size());
    combine_region_id=d1(rank_0.size());
    // archive_region_id=d1(rank_0.size());
    for(int i=0;i<rank_0.size();i++){
        archive[i]=tmp_archive[rank_0[i]];
        combine_region_id[i]=tmp_region_id[rank_0[i]];
        // archive_region_id[i]=tmp_region_id[rank_0[i]];
    }

}

void sede::create_rotMat(){
    dd2 temp_ma1(dim,dd1(dim,0));
    dd2 temp_ma2(dim,dd1(dim,0));
    for(int i=0;i<dim;i++)
        temp_ma1[i][i]=1;
    
    temp_ma1[dim-2][dim-1]=-sin(angle);
    temp_ma1[dim-2][dim-2]=cos(angle);
    temp_ma1[dim-1][dim-1]=cos(angle);
    temp_ma1[dim-1][dim-2]=sin(angle);
    
    for(int i=1;i<dim-1;i++)
        for(int j=0;j<=i;j++){
            for(int k=0;k<dim;k++){
                for(int l=0;l<dim;l++){
                    if(l!=(dim-1)-i-1 && l!=(dim-1)-j)
                        temp_ma2[k][l]=temp_ma1[l][l];
                    else
                        if(l==(dim-1)-i-1)
                            temp_ma2[k][l]=temp_ma1[(dim-1)-i-1][l]*cos(angle)+temp_ma1[(dim-1)-j][l]*sin(angle);
                        else if(l==(dim-1)-j)
                            temp_ma2[k][l]=temp_ma1[(dim-1)-i-1][l]*(-sin(angle))+temp_ma1[(dim-1)-j][l]*cos(angle);
                }
            }
            temp_ma1=temp_ma2;
        }
    
    rotMat=temp_ma1;
    rotMat_I=rotMat;
    for(int i=0;i<dim;i++)
        rotMat_I[i][i]=rotMat[i][i]-1;
}

void sede::spiral_tran(){
    double new_sol1=0.0;
    double new_sol2=0.0;
    dd3 temp_sample=sample;
    for(int i=0;i<searcher_num;i++)
        for (int j = 0; j < region_num; j++)
            for (int k = 0; k < sample_num; k++){
                // for(int ro_time=0;ro_time<10;ro_time++){
                    for(int l=0;l<dim;l++){
                        new_sol1=0.0;
                        new_sol2=0.0;
                        for(int m=0;m<dim;m++){
                            new_sol1+=rotMat[l][m]*sample[j][k][m];
                            new_sol2+=rotMat_I[l][m]*searcher[i][m];
                        }
                        sampleV[i][j][k][l]=new_sol1-new_sol2;
                        // temp_sample[j][k][l]=new_sol1-new_sol2;
                    }
                    // if(ro_time==9)
                    //     sampleV[i][j][k]=temp_sample[j][k];
                // }
                for(int l=0;l<dim;l++){
                    if(l<clip_bit_num){
                            if(sampleV[i][j][k][l]<region_bound[j][0][l]){
                                // sampleV[i][j][k][l]=(region_bound[j][0][l]+searcher[i][l])/2.0;
                                sampleV[i][j][k][l]=(region_bound[j][0][l]+sample[j][k][l])/2.0;
                            }
                            else if(sampleV[i][j][k][l]>region_bound[j][1][l]){
                                // sampleV[i][j][k][l]=(region_bound[j][1][l]+searcher[i][l])/2.0;
                                sampleV[i][j][k][l]=(region_bound[j][1][l]+sample[j][k][l])/2.0;
                            }
                        }
                        else{
                            if(sampleV[i][j][k][l]<lowerbound[l])
                                sampleV[i][j][k][l]=(lowerbound[l]+searcher[i][l])/2.0;
                            else if(sampleV[i][j][k][l]>upperbound[l])
                                sampleV[i][j][k][l]=(upperbound[l]+searcher[i][l])/2.0;
                    } 
                }
            }
    
}
void sede::de_mut(){
    double rnd;
    int num_player=sample_num*0.5;

    // if(tmp_iter>numIter*0.8)
    //     F=0.2;
    
    for (int i = 0; i < searcher_num; i++){
        for (int j = 0; j < region_num; j++){
            if(j==region_selected_by_s[i]){

                for (int k = 0; k < sample_num; k++){
                    //choose gbest1 and gbest2 with tournament selection
                    // int gbest1=rand()%sample_num;
                    int gbest2=rand()%sample_num;
                    int rand_s=rand()%sample_num;
                    while (searcher[i] == sample[j][gbest2])
                        gbest2 = rand()%sample_num;
                    // touranment selection
                    for (int l = 0; l < num_player - 1; l++)
                    {  
                        int r = rand()%sample_num;
                        if (sample_cdr_fit[j][r] < sample_cdr_fit[j][gbest2] && sample[j][r] != searcher[i])
                            gbest2 = r;
                    }
                    
                    
                    while (rand_s == gbest2 || sample[j][rand_s] == searcher[i])
                        rand_s = rand()%sample_num;
                    int d_rand=rand()%dim;
                    for(int l=0;l<dim;l++){
                        rnd=(double)rand()/RAND_MAX;
                        if(rnd<CR_i[j][k] ||d_rand==l){
                            sampleV[i][j][k][l]=sample[j][k][l]+F_i[j][k]*(searcher[i][l]-sample[j][k][l])
                                                +F_i[j][k]*(sample[j][gbest2][l]-sample[j][rand_s][l]);
    
                            //Keep the solution in the bound
                            if(l<clip_bit_num){
                                if(sampleV[i][j][k][l]<region_bound[j][0][l])
                                    sampleV[i][j][k][l]=(region_bound[j][0][l]+searcher[i][l])/2.0;
                                else if(sampleV[i][j][k][l]>region_bound[j][1][l])
                                    sampleV[i][j][k][l]=(region_bound[j][1][l]+searcher[i][l])/2.0;
                            }
                            else{
                                if(sampleV[i][j][k][l]<lowerbound[l])
                                    sampleV[i][j][k][l]=(lowerbound[l]+searcher[i][l])/2.0;
                                else if(sampleV[i][j][k][l]>upperbound[l])
                                    sampleV[i][j][k][l]=(upperbound[l]+searcher[i][l])/2.0;
                            }  
                                
                        }
                        else
                        {
                            sampleV[i][j][k][l]=sample[j][k][l];
                            // cout<<sampleV[i][j][k][l]<<" ";
                        }
                    }
                    // cout<<endl;
                    
                }
            }         
            else{
                for (int k = 0; k < sample_num; k++){
                    //choose gbest1 and gbest2 with tournament selection
                    // int gbest1=rand()%sample_num;
                    int gbest2=rand()%sample_num;
                    int rand_s1=rand()%sample_num;
                    int rand_s2=rand()%sample_num;
                    while (searcher[i] == sample[j][gbest2])
                        gbest2 = rand()%sample_num;
                    // touranment selection
                    for (int l = 0; l < num_player - 1; l++)
                    {
                        int r = rand()%sample_num;
                        // if (sample_cdr_fit[j][r] < sample_cdr_fit[j][gbest1] && r != rand_s1)
                        //     gbest1 = r;
                        if (sample_cdr_fit[j][r] < sample_cdr_fit[j][gbest2] && r != rand_s1 && searcher[i]!=sample[j][r])
                            gbest2 = r;
                    }
                    while (rand_s2 == gbest2 || rand_s2 == rand_s1)
                        rand_s2 = rand()%sample_num;

                    
                    int d_rand=rand()%dim;
                    for(int l=0;l<dim;l++){
                        rnd=(double)rand()/RAND_MAX;
                        if(rnd<CR_i[j][k] || d_rand==l){
                            sampleV[i][j][k][l]=sample[j][rand_s1][l]+F_i[j][k]*(searcher[i][l]-sample[j][rand_s1][l])
                                                +F_i[j][k]*(sample[j][gbest2][l]-sample[j][rand_s2][l]);
                            //Keep the solution in the bound
                            if(l<clip_bit_num){
                                if(sampleV[i][j][k][l]<region_bound[j][0][l]){
                                    // sampleV[i][j][k][l]=(region_bound[j][0][l]+searcher[i][l])/2.0;
                                    sampleV[i][j][k][l]=(region_bound[j][0][l]+sample[j][k][l])/2.0;
                                }
                                else if(sampleV[i][j][k][l]>region_bound[j][1][l]){
                                    // sampleV[i][j][k][l]=(region_bound[j][1][l]+searcher[i][l])/2.0;
                                    sampleV[i][j][k][l]=(region_bound[j][1][l]+sample[j][k][l])/2.0;
                                }
                            }
                            else{
                                if(sampleV[i][j][k][l]<lowerbound[l])
                                    sampleV[i][j][k][l]=(lowerbound[l]+searcher[i][l])/2.0;
                                else if(sampleV[i][j][k][l]>upperbound[l])
                                    sampleV[i][j][k][l]=(upperbound[l]+searcher[i][l])/2.0;
                            }  
                        }
                        else{
                            sampleV[i][j][k][l]=sample[j][k][l];
                        }
                    }
                    
                }
            }
        }
    }
        
}

//考慮三種計算CDR的方法
//1.將所有sample和sampleV加入archive並計算crowding_distance
/*2.每個點依次放入archive判斷單點在archive內的優異度，
以加入後的新擁擠距離與最小擁擠距離相減並加上新點之分佈平均*/
//3.利用Spread()公式判斷各點加入archive後的優異度
void sede::cal_cdr_fit(){

    combine_sol=archive;
    int archive_size=archive.size();
    d1 store_sample_pos;
    int sample_end=archive_size-1;
    
    //暫時將Sample加入archive內比較其優異度，若已有相同解在內則無須加入，紀錄位置即可
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            bool flag=true;
            // for(int k=0;k<archive_size;k++)
            //     if(combine_sol[k]==sample[i][j]){
            //         store_sample_pos.push_back(k);
            //         combine_region_id[k]=i;
            //         flag=false;
            //         break;
            //     }
            if(flag){
                combine_sol.push_back(sample[i][j]);
                sample_end++;
                store_sample_pos.push_back(sample_end);
                combine_region_id.push_back(i);
            }
        }
    
    //紀錄所有combine_sol的region_id
    /*
    combine_region_id=d1(combine_sol.size(),-1);
    
    for(int i=0;i<store_sample_pos.size();i++){
        combine_region_id[store_sample_pos[i]]=i/sample_num;
    }*/
    /*
    for(int i=0;i<combine_sol.size();i++){
        if(combine_region_id[i]==-1)
            for(int j=0;j<region_num;j++){
                if(equal(combine_sol[i][0].begin(),combine_sol[i][0].begin()+clip_bit_num,region_bit[j].begin())){
                    combine_region_id[i]=j;
                    break;
                }
            }
        
    }*/
    
    // cout<<"combine_sol_size1:"<<combine_sol.size()<<endl;
    //將sampleV暫時加入archive內比較他的優異度
    int sampleV_start=sample_end+1;
    // cout<<"SampleV_start:"<<sampleV_start<<endl;
    for(int i=0;i<region_num;i++)
        for(int j=0;j<searcher_num;j++)
            for(int k=0;k<sample_num;k++){
                combine_sol.push_back(sampleV[j][i][k]);
                // combine_region_id[sampleV_start]=i;
                combine_region_id.push_back(i);
                // combine_sol.push_back(sampleV[j][i][2*k+1]);
                // combine_region_id[sampleV_start+1]=i;
                // combine_region_id.push_back(i);
            }
    sampleV_start=sample_end+1;
    

    //取的所有二元解的實數解並取得Objective Value
    int combine_size=combine_sol.size();
    // dd2 combine_real_sol(combine_size);
    combine_fit.clear();//可能不需要
    // for(int i=0;i<combine_size;i++)
    //     combine_real_sol[i]=decode(combine_sol[i]);
    combine_fit=fitness(combine_sol);

    //紀錄dominated數量和dominate set
    // d2 set(combine_size,d1(1,0));
    // for(int i=0;i<combine_size-1;i++) 
    //     for(int j=i+1;j<combine_size;j++)
    //         if(combine_fit[i][0]<=combine_fit[j][0]&&combine_fit[i][1]<=combine_fit[j][1])
    //             if(combine_fit[i][0]==combine_fit[j][0]&&combine_fit[i][1]==combine_fit[j][1])
    //                 continue;
    //             else{
    //                 set[i].push_back(j);
    //                 set[j][0]++;
    //             }
    //         else if(combine_fit[i][0]>=combine_fit[j][0]&&combine_fit[i][1]>=combine_fit[j][1]){
    //             set[i][0]++;
    //             set[j].push_back(i);
    //         }
    //higher-objective
    d2 set(combine_size,d1(1,0));
    for(int i=0;i<combine_size-1;i++) 
        for(int j=i+1;j<combine_size;j++){
            int do_count=0;
            int equal_count=0;
            for(int ob=0;ob<obj_num;ob++)
                if(combine_fit[i][ob]<=combine_fit[j][ob]){
                    do_count++;
                    if(combine_fit[i][ob]==combine_fit[j][ob])
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
        for(int i=0;i<combine_size;i++)       
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
    
    dd2 use_for_cdr(rank.size());
    //做crowding_distance
    for(int i=0;i<rank.size();i++){
        if(rank[i].size()>2){
            // crowding_dis_assign(rank[i]);
            // crowding_entropy(rank[i]);
            use_for_cdr[i]=return_crowd_val(rank[i]);
        }
        else if(rank[i].size()==2){
            use_for_cdr[i].push_back(0.5);
            use_for_cdr[i].push_back(0.5);
        }
        else if(rank[i].size()==1)
            use_for_cdr[i].push_back(0.9);
        
    }
    
    //使用歐基里得距離
    // for(int i=0;i<rank.size();i++)
    //     if(rank[i].size()>2)
    //         Euclidean_crowding(rank[i]);
    
    // 1.嘗試填滿archive(bugfix:2020/03/18)
    if(rank[0].size()>max_pareto){
        
        rank_0=rank[0];
        //-----------------------------//
        // int delete_pos=rank[0].size()-max_pareto;
        // d1 delete_rank(delete_pos);
        // while(rank_0.size()>max_pareto){  
        //     crowding_dis_assign(rank_0);
        //     // crowding_entropy(rank_0);
        //     delete_rank[delete_pos-1]=rank_0[rank[0].size()-1];
        //     rank_0.pop_back();
        //     delete_pos--;
        // }
        
        // archive=dd2(rank_0.size());
        // d1 tmp_region_id=combine_region_id;
        // // archive_region_id=d1(rank_0.size());
        // for(int i=0;i<rank_0.size();i++){
        //     archive[i]=combine_sol[rank_0[i]];
        //     // archive_region_id[i]=tmp_region_id[rank_0[i]];
        // }
        // rank_0.insert(rank_0.end(),delete_rank.begin(),delete_rank.end());
        // rank[0]=rank_0;
        //-----------------------------//
        }
        else{
        rank_0=rank[0];
        int insert_rank_i=1;
        while(rank_0.size()<max_pareto && insert_rank_i<rank.size()){
            rank_0.insert(rank_0.end(),rank[insert_rank_i].begin(),rank[insert_rank_i].end());
            insert_rank_i++;
        }
        while(rank_0.size()>max_pareto)
            rank_0.pop_back();
        //-------------------------------//
        // archive=dd2(rank_0.size());
        // d1 tmp_region_id=combine_region_id;
        // combine_region_id=d1(rank_0.size());
        // // archive_region_id=d1(rank_0.size());
        // for(int i=0;i<rank_0.size();i++){
        //     archive[i]=combine_sol[rank_0[i]];
        //     combine_region_id[i]=tmp_region_id[rank_0[i]];
        //     // archive_region_id[i]=tmp_region_id[rank_0[i]];
        // }

        //-------------------------------//
        }
    
    //2.只存rank1的解但可能會出錯
    // rank_0=rank[0];
    
    ///----------------------------------------------------///
    
    combine_region_id_2d=d2(region_num);
    for(int i=0;i<rank.size();i++)
        for(int j=0;j<rank[i].size();j++)
            combine_region_id_2d[combine_region_id[rank[i][j]]].push_back(rank[i][j]);
    
    //計算所有解的CDR
    // double giving_level=(double)combine_size;
    // double rank_size=(double)rank.size();
    // combine_cdr_fit.assign(combine_size,0.0);
    // for(int i=0;i<rank_size;i++)
    //     for(int j=0;j<rank[i].size();j++){
    //         combine_cdr_fit[rank[i][j]]=giving_level*pow(dominate_pressure,(rank_size-i));
    //         giving_level=giving_level-1;
    //     }

    //將CDR值正規化
    // double giving_level=(double)combine_size;
    // double rank_size=(double)rank.size();
    // combine_cdr_fit.assign(combine_size,0.0);
    // for(int i=0;i<rank_size;i++)
    //     for(int j=0;j<rank[i].size();j++){
    //         combine_cdr_fit[rank[i][j]]=((giving_level-1)/(combine_size-1))*pow(dominate_pressure,((rank_size-i)-1)/(rank_size-1));
    //         giving_level=giving_level-1;
    //     }
    
    //針對crowding distance計算cdr
    double rank_size=(double)rank.size();
    combine_cdr_fit.assign(combine_size,0.0);
    double d_p=20.0;
    for(int i=0;i<rank_size;i++)
        for(int j=0;j<rank[i].size();j++)
            combine_cdr_fit[rank[i][j]]=use_for_cdr[i][j]+d_p*(rank_size-i);
            
    
    //存入Sample_cdr_fit和SampleV_cdr_fit
    int count=0;
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            sample_cdr_fit[i][j]=combine_cdr_fit[store_sample_pos[count]];
            count++;
        }
    count=sampleV_start;

    for(int i=0;i<region_num;i++)
        for(int j=0;j<searcher_num;j++)
            for(int k=0;k<sample_num;k++){
                sampleV_cdr_fit[j][i][k]=combine_cdr_fit[count];
                count++;
            }
    
}
sede::dd1 sede::return_crowd_val(d1& this_rank){
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

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
        for (int i=0; i<obj_num; i++){
            init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
            init_dis[sorting_id[i][front_size-1]]= std::numeric_limits<double>::infinity();

        }
        for(int i=0;i<obj_num;i++)
            for(int j=1;j<front_size-1;j++){
                // if(init_dis[j]!=numeric_limits<double>::infinity()){
                    // if(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]== combine_fit[this_rank[sorting_id[i][0]]][i])
                    //     init_dis[sorting_id[i][j]]+=0.0;
                    // else
                        init_dis[sorting_id[i][j]]+=(combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i])/(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
                // }
            }
            //距離密度值進行正規化
            for (int j=0; j<front_size; j++){
                if (init_dis[j]!= numeric_limits<double>::infinity())
                    init_dis[j] = init_dis[j]/obj_num;
            }

        // bool check_inf=true;
        // if(obj_num>=3)
        //     for(int i=0;i<front_size;i++)
        //         if(init_dis[i]!=inf)
        //             check_inf=false;
        // if(!check_inf){
            //two_times infinity
            double max_value=-numeric_limits<double>::infinity();
            for(int i=0;i<front_size;i++)
                if(init_dis[i]>max_value && init_dis[i]!=numeric_limits<double>::infinity())
                    max_value=init_dis[i];
            for(int i=0;i<front_size;i++)
                if(init_dis[i]==numeric_limits<double>::infinity())
                    init_dis[i]=max_value*2;
        // }
        // else
        //     for(int i=0;i<front_size;i++)
        //         init_dis[i]=1.0/front_size;
    }
    //when front size = 2
    else if(front_size==2){
        init_dis[0]=0.5;
        init_dis[1]=0.5;
    }
    else if(front_size==1)
        init_dis[0]=0.9;
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end());
    reverse(init_dis.begin(),init_dis.end());
    return init_dis;
}
//Generate new crossover_rate and F value(Need to fine tune)
void sede::generate_new_Cr_F(){
    random_device rd;
    mt19937 generator(rd());
    normal_distribution<double> cr_dist(u_cr,0.1);
    normal_distribution<double> f_dist(u_f,0.1);

    default_random_engine gen = std::default_random_engine(rd());
    uniform_real_distribution<double> rndF(0,1.2);
    
    int one_third=(region_num*sample_num)/3;
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){

            CR_i[i][j]=cr_dist(generator);
            F_i[i][j]=f_dist(generator);
            // if(((double)rand()/RAND_MAX)<0.1)
            //     CR_i[i][j]=0.0+((double)rand()/RAND_MAX)*(0.9-0.0);
            
            // if(((double)rand()/RAND_MAX)<0.1)
            //     F_i[i][j]=0.1+((double)rand()/RAND_MAX)*(0.9-0.1);
        }
    for(int i=0;i<region_num;i++){
        if(one_third!=0){
            F_i[i][rand()%sample_num]=rndF(gen);
            one_third--;
        }
        else
            break;
        if(i==region_num-1)
            i=(i+1)%region_num;   
    }
    S_cr.clear();
    S_f.clear();
}
void sede::update_uCr_uF(){
    double c=0.9;
    //
    double mean_S_cr=0;
    for(int i=0;i<S_cr.size();i++)
        mean_S_cr+=S_cr[i];
    if(mean_S_cr!=0)
        mean_S_cr/=S_cr.size();
    else
        mean_S_cr=0;
    u_cr=(1-c)*u_cr+c*mean_S_cr;
    //
    double lemer=0.0;
    double sum1=0,sum2=0;
    for(int i=0;i<S_f.size();i++){
        sum1+=pow(S_f[i],2.0);
        sum2+=S_f[i];
    }
    if(sum2!=0)
        lemer=sum1/sum2;
    else
        lemer=0; 
    
    u_f=(1-c)*u_f+c*lemer;
}
void sede::combine_region(){
    clip_bit_num--;
    region_num/=2;
    divide_region();
    region_num*=2;
    // sample_num*=2;
    
    dd3 new_sample(region_num/2,dd2(sample_num*2,dd1(dim,0.0)));
    dd2 new_sample_fit(region_num/2,dd1(sample_num*2));
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            new_sample[i/2][j+(i%2)*sample_num]=sample[i][j];
            new_sample_fit[i/2][j+(i%2)*sample_num]=sample_cdr_fit[i][j];
        }
    
    //(1)Reduce the number of Samples and increase the number of Searcher(This method may can increase convergence)//
    // d2 sort_id(region_num/2,d1(sample_num*2));
    // dd3 half_best_sample(region_num/2,dd2(sample_num));
    // dd2 half_best_sample_fit(region_num/2,dd1(sample_num));
    // for(int i=0;i<region_num/2;i++)
    //     for(int j=0;j<sample_num*2;j++)
    //         sort_id[i][j]=j;
    // for(int i=0;i<region_num/2;i++)
    //     quick_sort(new_sample_fit[i],sort_id[i],0,sample_num*2-1);

    // for(int i=0;i<region_num/2;i++)
    //     for(int j=0;j<sample_num;j++){
    //         half_best_sample[i][j]=new_sample[i][sort_id[i][j]];
    //         half_best_sample_fit[i][j]=new_sample_fit[i][j];
    //     }
        
    //-------------------------------------------------------------------------------------------------------//
    dd4 new_sampleV(searcher_num,dd3(region_num/2,dd2(sample_num*2,dd1(dim,0.0))));
    dd3 new_sampleV_fit(searcher_num,dd2(region_num/2,dd1(sample_num*2)));
    for(int i=0;i<searcher_num;i++)
        for(int j=0;j<region_num;j++){
            for(int k=0;k<sample_num;k++){
                new_sampleV[i][j/2][k+(j%2)*sample_num]=sampleV[i][j][k];
                new_sampleV_fit[i][j/2][k+(j%2)*sample_num]=sampleV_cdr_fit[i][j][k];
            }
        }
    
    dd2 new_F_i(region_num/2,dd1(sample_num*2));
    dd2 new_Cr_i(region_num/2,dd1(sample_num*2));
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            new_F_i[i/2][j+(i%2)*sample_num]=F_i[i][j];
            new_Cr_i[i/2][j+(i%2)*sample_num]=CR_i[i][j];
        }
   
    //searcher_select_region
    for(int i=0;i<searcher_num;i++)
        region_selected_by_s[i]=region_selected_by_s[i]/2;
     
    //---------------------------------(1)Increase Searcher------------------------------------//

    // dd2 new_searcher=searcher;
    // new_searcher.insert(new_searcher.end(),searcher.begin(),searcher.end());
    // searcher=new_searcher;
    // d1 new_belong=region_selected_by_s;
    // new_belong.insert(new_belong.end(),region_selected_by_s.begin(),region_selected_by_s.end());
    // region_selected_by_s=new_belong;

    // dd2 half_best_F(region_num/2,dd1(sample_num));
    // dd2 half_best_Cr(region_num/2,dd1(sample_num));
    // for(int i=0;i<region_num/2;i++)
    //     for(int j=0;j<sample_num;j++){
    //         new_F_i[i][j]=new_F_i[i][sort_id[i][j]];
    //         new_Cr_i[i][j]=new_Cr_i[i][sort_id[i][j]];
    //     }
    
    //--------------------------------------------------------------------------------------//
    //id
    for(int i=0;i<combine_region_id.size();i++){
        combine_region_id[i]=combine_region_id[i]/2;
        // archive_region_id[i]=archive_region_id[i]/2;
    }
    //ta and tb init
    tb.assign(region_num,0.0);
    ta.assign(region_num,1.0);
    //sample and sampleV init
    //----------------------------//
    sample=new_sample;
    sample_cdr_fit=new_sample_fit;
    // sample=half_best_sample;
    // sample_cdr_fit=half_best_sample_fit;
    //----------------------------//
    F_i=new_F_i;
    CR_i=new_Cr_i;
    //-------(1)---------//
    sampleV=new_sampleV;
    sampleV_cdr_fit=new_sampleV_fit;
    sample_num*=2;
    // sampleV.assign(searcher_num*2,dd3(region_num/2,dd2(sample_num,dd1(dim,0.0))));
    // sampleV_cdr_fit.assign(searcher_num*2,dd2(region_num/2,dd1(sample_num)));
    // searcher_num=searcher_num*2;
    //------------------//
    region_num/=2;
    
}
//--------------------------------------------------------------------//(試試是否正確)
void sede::crowding_dis_assign(d1& this_rank){
   
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    
    for(int i=0;i<obj_num;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    }
    
    //---------2目標專用---------//
    
    // for(int j=0;j<front_size;j++)
    //     sorting_id[0][j]=j;
    // obj_id[0]=this_rank;
    // quick_sort_obj(0,obj_id[0],sorting_id[0],0,front_size-1);
    // obj_id[1]=obj_id[0];
    // sorting_id[1]=sorting_id[0];
    // reverse(obj_id[1].begin(),obj_id[1].end());
    // reverse(sorting_id[1].begin(),sorting_id[1].end());
    
    //------------------------//

    //計算擁擠程度
    for (int i=0; i<obj_num; i++){
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
        init_dis[sorting_id[i][front_size-1]]= std::numeric_limits<double>::infinity();
    }
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<front_size-1;j++){
            // if(init_dis[j]!=numeric_limits<double>::infinity()){
                // if(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]== combine_fit[this_rank[sorting_id[i][0]]][i])
                //     init_dis[sorting_id[i][j]]+=0.0;
                // else
                    init_dis[sorting_id[i][j]]+=(combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i])/(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
            // }
        }
    //距離密度值進行正規化
    for (int j=0; j<front_size; j++){
        if (init_dis[j]!= numeric_limits<double>::infinity())
           init_dis[j] = init_dis[j]/obj_num;
    }
    //
    // d1 tmp_rank=this_rank;
    // quick_sort(init_dis,tmp_rank,0,front_size-1);
    // deleted_pos=tmp_rank[0];
    
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end()); 
}
void sede::crowding_entropy(d1& this_rank){
   int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    
    for(int i=0;i<obj_num;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    }
    //計算擁擠程度
    for (int i=0; i<obj_num; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
    double E_ij,c_ij,pl_ij,pu_ij;
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<front_size-1;j++){
            if(init_dis[sorting_id[i][j]]!=numeric_limits<double>::infinity()){
                
                c_ij=combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i];
                
                pu_ij=(combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j]]][i])/c_ij;
                
                pl_ij=(combine_fit[this_rank[sorting_id[i][j]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i])/c_ij;
                E_ij=-(pl_ij*log2(pl_ij)+pu_ij*log2(pu_ij));
                // cout<<E_ij<<endl;;
                if(c_ij==0||pu_ij==0||pl_ij==0)
                    init_dis[sorting_id[i][j]]+=0.0;
                else
                    init_dis[sorting_id[i][j]]+=(c_ij*E_ij)/(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
                
            }    
        }
    for(int i=0;i<init_dis.size();i++)
        if(init_dis[i]!=std::numeric_limits<double>::infinity())
            init_dis[i]=init_dis[i]*(-1);
    
    // for(int i=0;i<init_dis.size();i++)
    //     cout<<init_dis[i]<<endl;
    // cout<<endl;
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end()); 
    
}
void sede::cal_ev(){

    //暫存所有archive,sample,sampleV並做fast non-dominated sort記錄所有的rank和crowding distance
    //
    
    cal_cdr_fit();
    
    double Tj,Vj,Mj;

    double sample_fit_sum=0.0;
    //計算M_j(考慮兩種作法:1.取目前sample最好之值 2.取archive內相同區域且最好的值)
    //以下使用方法為(1)
    for(int i=0;i<region_num;i++){
        double tmp_best=-1;
        int best_pos=-1;
        for(int j=0;j<sample_num;j++){
            sample_fit_sum+=sample_cdr_fit[i][j];
            if(sample_cdr_fit[i][j]>tmp_best){
                tmp_best=sample_cdr_fit[i][j];
                best_pos=j;
            }
            region_best_fit[i]=tmp_best;
        }
    }

    //計算V_ij
    dd2 Vij(searcher_num,dd1(region_num,0.0));
     for (int i = 0; i < searcher_num; i++) {
        for (int j = 0; j < region_num; j++) {
            for (int k = 0; k < sample_num; k++) {
                Vij[i][j] += sampleV_cdr_fit[i][j][k];
            }
            Vij[i][j] /= sample_num; 
        }
    }
     EV.assign(searcher_num,dd1(region_num,0.0));
    //計算所有Searcher各區域的期望值
     for (int i = 0; i < searcher_num; i++) 
        for (int j = 0; j < region_num; j++) 
            EV[i][j] = (tb[j]/ta[j])*Vij[i][j]*(region_best_fit[j]/sample_fit_sum);
   
    // for(int i=0;i<searcher_num;i++)
    //     searcher_fit[i]=evaluation(searcher[i]);
    choose_sample(1);
}

void sede::choose_sample(int choose_method){
    //考慮三種作法:
    //1.sample與對應sampleV比較並取最大值
    //2.同區域內取前幾大(sample數)做比較並替換原sample
    //3.採用Tournament Selection來選Sample
    //4.clonal selection
    //以下採用(1)
    if(choose_method==1){
        for (int i = 0; i < searcher_num ; i++) 
            for(int j = 0; j < region_num; j++) 
                for (int k = 0; k < sample_num; k++) {
                    if (sampleV_cdr_fit[i][j][k] > sample_cdr_fit[j][k]) {
                        // for (int l = clip_bit_num; l < numPatterns; l++)
                        sample[j][k]= sampleV[i][j][k];
                        sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k];
                        //
                        S_cr.push_back(CR_i[j][k]);
                        S_f.push_back(F_i[j][k]);
                    }
                    // if (sampleV_cdr_fit[i][j][k*2+1] > sample_cdr_fit[j][k]) {
                    //     // for (int l = clip_bit_num; l < numPatterns; l++)
                    //     sample[j][k] = sampleV[i][j][k*2+1];
                    //     sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k*2+1];
                    // }
                }
    }
    //以下採用(2)(未寫)
    else if(choose_method==2){

    }
    //以下採用(3)(2020/03/17 add in)
    else if(choose_method==3){
        int selected_samV1,selected_samV2,V1_1,V1_3,V2_1,V2_3;
        for(int i=0;i<region_num;i++){
            for(int j=0;j<sample_num;j++){
                selected_samV1=rand()%(sample_num*2*searcher_num);
                do{
                    selected_samV2=rand()%(sample_num*2*searcher_num);
                }while(selected_samV2==selected_samV1);
                V1_1=selected_samV1/(sample_num*2);
                V1_3=selected_samV1%(sample_num*2);
                V2_1=selected_samV2/(sample_num*2);
                V2_3=selected_samV2%(sample_num*2);
                // cout<<V1_1<<" "<<V1_3<<endl;
                // cout<<V2_1<<" "<<V2_3<<endl<<endl;
               
                if(sampleV_cdr_fit[V1_1][i][V1_3]>sampleV_cdr_fit[V2_1][i][V2_3]){
                    //以下if考慮是否刪除增加多樣性
                    if(sample_cdr_fit[i][j]<sampleV_cdr_fit[V1_1][i][V1_3]){
                        sample[i][j]= sampleV[V1_1][i][V1_3];
                        sample_cdr_fit[i][j] = sampleV_cdr_fit[V1_1][i][V1_3];
                    }
                }
                else{
                    //以下if考慮是否刪除增加多樣性
                    if(sample_cdr_fit[i][j]<sampleV_cdr_fit[V2_1][i][V2_3]){
                        sample[i][j]= sampleV[V2_1][i][V2_3];
                        sample_cdr_fit[i][j] = sampleV_cdr_fit[V2_1][i][V2_3];
                    }
                }
                
            }
        }
    }
    else if(choose_method==4){

    }
}
void sede::select_player(){

    for (int i = 0; i < region_num; i++)
        tb[i]++;
    
    int searcher_count=0;
    d1 s_candidate(searcher_num);
    dd1 s_candidate_val(searcher_num);
    d1 s_candidate_region(searcher_num);

    d1 region_best_sol(region_num);
    for(int i=0;i<region_num;i++)
        region_best_sol[i]=combine_region_id_2d[i][0];

    for (int i = 0; i < searcher_num; i++) {
        
        int first_select = rand() % region_num;
        double first_ev = EV[i][first_select];
        // if(tmp_iter<numIter*0.5)
            for(int j=0;j<tour_size;j++){
                int second_select= rand() % region_num;
                if (EV[i][second_select] >  first_select) {
                    first_select = second_select;
                    first_ev = EV[i][first_select];
                }
            }
        // else{
        //     first_select=0;
        //     for(int j=1;j<region_num;j++)
        //         if(EV[i][j]>EV[i][first_select])
        //             first_select=j;
        // }
        // for(int j=0;j<region_num;j++)
        //     if(EV[i][j]>EV[i][first_select]){
        //         first_ev=EV[i][j];
        //         first_select=j;
        //     }
        //決定每個Searcher，有兩種作法
        //1.選取最好------->(1)選前幾好的(X)，(2)一律選最好的(V)
        //2.使用Tournament Selection選取
        //以下採用(1.1)
        if(combine_region_id_2d[first_select].size()>0){
            //--------------Use for clonal selection----------------//
            // s_candidate[searcher_count]=combine_region_id_2d[first_select][0];
            // s_candidate_val[searcher_count]=combine_cdr_fit[combine_region_id_2d[first_select][0]];
            // s_candidate_region[searcher_count]=first_select;
            // searcher_count++;
            //-----------------Original-----------------//
            searcher[i]=combine_sol[combine_region_id_2d[first_select][0]];
            //將以下這行註解後即為(1.2)的方法
            combine_region_id_2d[first_select].erase(combine_region_id_2d[first_select].begin());
        }
        else{
            double tmp_best=-1;
            int best_pos=-1;
            for (int j = 0; j < sample_num; j++) 
                if (sample_cdr_fit[first_select][j] >tmp_best) {
                    tmp_best=sample_cdr_fit[first_select][j];
                    best_pos=j;
                }
            //------------Original------------//
            searcher[i]=sample[first_select][best_pos];

            //-----------Use for clonal selection-------------//
            // s_candidate[searcher_count]=region_best_sol[first_select];
            // s_candidate_val[searcher_count]=combine_cdr_fit[region_best_sol[first_select]];
            // s_candidate_region[searcher_count]=first_select;
            // searcher_count++;
            //--------select best in that region--------//
            // searcher[i]=combine_sol[region_best_sol[first_select]];
        }
        region_selected_by_s[i]=first_select;
        // update region_it[i] and region_hl[i];
        ta[first_select]++;
        tb[first_select] = 1;
    }

    //------------clonal selection--------------//
    // dd1 tmp_s_candidate_val=s_candidate_val;
    // quick_sort(s_candidate_val,s_candidate,0,searcher_num-1);
    // quick_sort(tmp_s_candidate_val,s_candidate_region,0,searcher_num-1);
    // int NC_num=5;
    // d1 good_set(NC_num);
    // dd1 good_set_val(NC_num);
    // d1 good_set_region(NC_num);
    // double good_sum=0;
    // for(int i=0;i<NC_num;i++){
    //     good_set[i]=s_candidate[searcher_num-1-i];
    //     good_set_val[i]=s_candidate_val[searcher_num-1-i];
    //     good_set_region[i]=s_candidate_region[searcher_num-1-i];
    //     good_sum+=s_candidate_val[searcher_num-1-i];
    // }
    // d1 new_searcher_set;
    // int set_num;
    // region_selected_by_s.clear();
    // for(int i=0;i<NC_num;i++){
    //     set_num=searcher_num*(good_set_val[i]/good_sum);
    //     for(int j=0;j<set_num;j++){
    //         new_searcher_set.push_back(good_set[i]);
    //         region_selected_by_s.push_back(good_set_region[i]);
    //     }
    // }
    // int fix_size=searcher_num-new_searcher_set.size();
    // int random_select;
    // for(int i=0;i<fix_size;i++){
    //     random_select=rand()%NC_num;
    //     new_searcher_set.push_back(good_set[random_select]);
    //     region_selected_by_s.push_back(good_set_region[random_select]);
    // }
    // for(int i=0;i<searcher_num;i++)
    //     searcher[i]=combine_sol[new_searcher_set[i]];
    //-----------------------------------------------//
}
void sede::archive_delete(){
        
        double r;
        r=(double)rand()/RAND_MAX;
        // if(tmp_iter<(double)numIter*0.85)
        //     decom_archive_delete();
        // else
            while(rank_0.size()>max_pareto){
                
                // if(r<0.7)
                    crowding_dis_assign(rank_0);
                    // crowding_entropy(rank_0);
                // else
                //     Euclidean_crowding(rank_0);
                
                rank_0.pop_back();
            }
            
        archive=dd2(rank_0.size());
        d1 tmp_region_id=combine_region_id;
        combine_region_id=d1(rank_0.size());
        for(int i=0;i<rank_0.size();i++){
            archive[i]=combine_sol[rank_0[i]];
            combine_region_id[i]=tmp_region_id[rank_0[i]];
        }

}
void sede::decom_archive_delete(){
    //以TCH值大小判斷欲刪除的點
    if(rank_0.size()>max_pareto){
        double min_tch,cur_tch;
        int weight_belong;
        d1 new_archive;
        d2 weight_belong_set(weight_vector_size);
        dd2 weight_belong_tch(weight_vector_size);
        for(int i=0;i<rank_0.size();i++){
            min_tch=inf;
            for(int j=0;j<weight_vector_size;j++){
                cur_tch=TCH(combine_fit[rank_0[i]],weight_vector[j]);
                if(cur_tch<min_tch){
                    min_tch=cur_tch;
                    weight_belong=j;
                }
            }
            weight_belong_set[weight_belong].push_back(rank_0[i]);
            weight_belong_tch[weight_belong].push_back(min_tch);
        }
        int min_pos;
        for(int i=0;i<weight_vector_size;i++){
            if(weight_belong_set[i].size()>1){
                min_tch=weight_belong_tch[i][0];
                min_pos=0;
                for(int j=1;j<weight_belong_set[i].size();j++){
                    if(weight_belong_tch[i][j]<min_tch){
                        min_tch=weight_belong_tch[i][j];
                        min_pos=j;
                    }
                }
                new_archive.push_back(weight_belong_set[i][min_pos]);
            }
        }
        rank_0=new_archive;
    }

}
double sede::TCH(dd1& sol_fit,dd1& cur_wei){
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
void sede::Euclidean_crowding(d1& this_rank){
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    // for(int i=0;i<obj_num;i++){
    //     for(int j=0;j<front_size;j++)
    //         sorting_id[i][j]=j;
    //     obj_id[i]=this_rank;
    //     quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    // }
    //---------2目標專用---------//
    for(int j=0;j<front_size;j++)
        sorting_id[0][j]=j;
    obj_id[0]=this_rank;
    quick_sort_obj(0,obj_id[0],sorting_id[0],0,front_size-1);
    obj_id[1]=obj_id[0];
    sorting_id[1]=sorting_id[0];
    reverse(obj_id[1].begin(),obj_id[1].end());
    reverse(sorting_id[1].begin(),sorting_id[1].end());
    //------------------------//

    //用歐基里德距離計算擁擠度
    init_dis[sorting_id[0][0]]=std::numeric_limits<double>::infinity();
    init_dis[sorting_id[0][front_size-1]]=std::numeric_limits<double>::infinity();

    for(int i=1;i<front_size-1;i++){
        init_dis[sorting_id[0][i]]=
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i-1]]])+
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i+1]]]);
    }
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end());
}
void sede::Euclidean_archive_d(d1& this_rank){
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    // for(int i=0;i<obj_num;i++){
    //     for(int j=0;j<front_size;j++)
    //         sorting_id[i][j]=j;
    //     obj_id[i]=this_rank;
    //     quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    // }
    //---------2目標專用---------//
    for(int j=0;j<front_size;j++)
        sorting_id[0][j]=j;
    obj_id[0]=this_rank;
    quick_sort_obj(0,obj_id[0],sorting_id[0],0,front_size-1);
    obj_id[1]=obj_id[0];
    sorting_id[1]=sorting_id[0];
    reverse(obj_id[1].begin(),obj_id[1].end());
    reverse(sorting_id[1].begin(),sorting_id[1].end());
    //------------------------//
    //計算擁擠程度
    init_dis[sorting_id[0][0]]=std::numeric_limits<double>::infinity();
    init_dis[sorting_id[0][front_size-1]]=std::numeric_limits<double>::infinity();

    for(int i=1;i<front_size-1;i++){
        init_dis[sorting_id[0][i]]=
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i-1]]])+
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i+1]]]);
    }

    double min_sol_dis;
    int min_sol_pos;
    d1 min_sorting_id(obj_num);
    
    while(front_size>max_pareto){
        min_sol_dis=numeric_limits<double>::infinity();
        min_sol_pos=-1;
        for(int i=0;i<front_size;i++)
            if(min_sol_dis>init_dis[i]){
                // cout<<"in+min"<<endl;
                min_sol_dis=init_dis[i];
                min_sol_pos=i;
            }
        
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]==min_sol_pos)
                {
                    // cout<<i<<":"<<j<<endl;
                    if(init_dis[sorting_id[i][j+1]]!=numeric_limits<double>::infinity())
                        init_dis[sorting_id[i][j+1]]=0.0;
                    if(init_dis[sorting_id[i][j-1]]!=numeric_limits<double>::infinity())  
                        init_dis[sorting_id[i][j-1]]=0.0;
                    min_sorting_id[i]=j;
                    break;
                }
        for(int i=0;i<obj_num;i++){
            
            if(init_dis[sorting_id[i][min_sorting_id[i]+1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]+1]]=
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+2]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]])+
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]);
                
            }
            if(init_dis[sorting_id[i][min_sorting_id[i]-1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]-1]]=
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-2]]])+
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]]);
            }
            // cout<<i<<"="<<init_dis[sorting_id[i][min_sorting_id[i]+1]]<<" "<<init_dis[sorting_id[i][min_sorting_id[i]-1]]<<endl;
        }
        // cout<<front_size<<endl;
        // cout<<endl;
        front_size--;
        for(int i=0;i<obj_num;i++)
            sorting_id[i].erase(sorting_id[i].begin()+min_sorting_id[i]);
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]>min_sol_pos)
                    sorting_id[i][j]--;
        this_rank.erase(this_rank.begin()+min_sol_pos);
        init_dis.erase(init_dis.begin()+min_sol_pos);
        
    }
    
    archive=dd2(rank_0.size());
    d1 tmp_region_id=combine_region_id;
    combine_region_id=d1(rank_0.size());
    for(int i=0;i<rank_0.size();i++){
        archive[i]=combine_sol[rank_0[i]];
        combine_region_id[i]=tmp_region_id[rank_0[i]];
    }
}
//計算cdr用
void sede::Euclidean_d_evaluate(d1& this_rank){
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);
    d1 deleted_sol;
    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    // for(int i=0;i<obj_num;i++){
    //     for(int j=0;j<front_size;j++)
    //         sorting_id[i][j]=j;
    //     obj_id[i]=this_rank;
    //     quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    // }
    //---------2目標專用---------//
    for(int j=0;j<front_size;j++)
        sorting_id[0][j]=j;
    obj_id[0]=this_rank;
    quick_sort_obj(0,obj_id[0],sorting_id[0],0,front_size-1);
    obj_id[1]=obj_id[0];
    sorting_id[1]=sorting_id[0];
    reverse(obj_id[1].begin(),obj_id[1].end());
    reverse(sorting_id[1].begin(),sorting_id[1].end());
    //------------------------//
    //計算擁擠程度
    init_dis[sorting_id[0][0]]=std::numeric_limits<double>::infinity();
    init_dis[sorting_id[0][front_size-1]]=std::numeric_limits<double>::infinity();

    for(int i=1;i<front_size-1;i++){
        init_dis[sorting_id[0][i]]=
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i-1]]])+
        Euclidean_dis(combine_fit[this_rank[sorting_id[0][i]]],combine_fit[this_rank[sorting_id[0][i+1]]]);
    }

    double min_sol_dis;
    int min_sol_pos;
    d1 min_sorting_id(obj_num);
    
    while(front_size>max_pareto){
        min_sol_dis=numeric_limits<double>::infinity();
        min_sol_pos=-1;
        for(int i=0;i<front_size;i++)
            if(min_sol_dis>init_dis[i]){
                
                min_sol_dis=init_dis[i];
                min_sol_pos=i;
            }
        deleted_sol.push_back(this_rank[min_sol_pos]);
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]==min_sol_pos)
                {
                    if(init_dis[sorting_id[i][j+1]]!=numeric_limits<double>::infinity())
                        init_dis[sorting_id[i][j+1]]=0.0;
                    if(init_dis[sorting_id[i][j-1]]!=numeric_limits<double>::infinity())  
                        init_dis[sorting_id[i][j-1]]=0.0;
                    min_sorting_id[i]=j;
                    break;
                }
        for(int i=0;i<obj_num;i++){
            
            if(init_dis[sorting_id[i][min_sorting_id[i]+1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]+1]]=
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+2]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]])+
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]);
                
            }
            if(init_dis[sorting_id[i][min_sorting_id[i]-1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]-1]]=
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-2]]])+
                Euclidean_dis(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]]
                ,combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]]);
            }
            // cout<<i<<"="<<init_dis[sorting_id[i][min_sorting_id[i]+1]]<<" "<<init_dis[sorting_id[i][min_sorting_id[i]-1]]<<endl;
        }
        // cout<<front_size<<endl;
        // cout<<endl;
        front_size--;
        for(int i=0;i<obj_num;i++)
            sorting_id[i].erase(sorting_id[i].begin()+min_sorting_id[i]);
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]>min_sol_pos)
                    sorting_id[i][j]--;
        this_rank.erase(this_rank.begin()+min_sol_pos);
        init_dis.erase(init_dis.begin()+min_sol_pos);
        
    }
    
    archive=dd2(rank_0.size());
    d1 tmp_region_id=combine_region_id;
    combine_region_id=d1(rank_0.size());
    for(int i=0;i<rank_0.size();i++){
        archive[i]=combine_sol[rank_0[i]];
        combine_region_id[i]=tmp_region_id[rank_0[i]];
    }
    
    int deleted_sol_size=deleted_sol.size();
    // cout<<deleted_sol_size<<endl;
   
    dd1 deleted_init_dis(deleted_sol_size,0.0);
    for(int i=0;i<deleted_sol_size;i++)
        for(int j=0;j<front_size;j++)
            if(combine_fit[deleted_sol[i]][0]<combine_fit[this_rank[sorting_id[0][j]]][0]){
                deleted_init_dis[i]=
                Euclidean_dis(combine_fit[this_rank[sorting_id[0][j]]]
                ,combine_fit[deleted_sol[i]])+
                Euclidean_dis(combine_fit[this_rank[sorting_id[0][j-1]]]
                ,combine_fit[deleted_sol[i]]);
                break;
            }
        
    quick_sort(deleted_init_dis,deleted_sol,0,deleted_sol_size-1);
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(deleted_sol.begin(),deleted_sol.end());
    reverse(this_rank.begin(),this_rank.end());

    this_rank.insert(this_rank.end(),deleted_sol.begin(),deleted_sol.end());

}
//計算cdr用
void sede::crowding_d_evaluate(d1& this_rank){
    
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);
    d1 deleted_sol;
    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    //--------各目標均可用--------//
    // for(int i=0;i<obj_num;i++){
    //     for(int j=0;j<front_size;j++)
    //         sorting_id[i][j]=j;
    //     obj_id[i]=this_rank;
    //     quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    // }
    
    //---------2目標專用---------//
    for(int j=0;j<front_size;j++)
        sorting_id[0][j]=j;
    obj_id[0]=this_rank;
    quick_sort_obj(0,obj_id[0],sorting_id[0],0,front_size-1);
    obj_id[1]=obj_id[0];
    sorting_id[1]=sorting_id[0];
    reverse(obj_id[1].begin(),obj_id[1].end());
    reverse(sorting_id[1].begin(),sorting_id[1].end());
    //------------------------//
    //計算擁擠程度
    for (int i=0; i<obj_num; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<front_size-1;j++){
            init_dis[sorting_id[i][j]]+=(combine_fit[this_rank[sorting_id[i][j+1]]][i]
            -combine_fit[this_rank[sorting_id[i][j-1]]][i])/
            (combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
        }

    double min_sol_dis;
    int min_sol_pos;
    d1 min_sorting_id(obj_num);
    
    while(front_size>max_pareto){
        min_sol_dis=numeric_limits<double>::infinity();
        min_sol_pos=-1;
        for(int i=0;i<front_size;i++)
            if(min_sol_dis>init_dis[i]){
                min_sol_dis=init_dis[i];
                min_sol_pos=i;
            }
        deleted_sol.push_back(this_rank[min_sol_pos]);
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]==min_sol_pos)
                {
                    // cout<<i<<":"<<j<<endl;
                    if(init_dis[sorting_id[i][j+1]]!=numeric_limits<double>::infinity())
                        init_dis[sorting_id[i][j+1]]=0.0;
                    if(init_dis[sorting_id[i][j-1]]!=numeric_limits<double>::infinity())  
                        init_dis[sorting_id[i][j-1]]=0.0;
                    min_sorting_id[i]=j;
                    break;
                }
                
        for(int i=0;i<obj_num;i++){
            
            if(init_dis[sorting_id[i][min_sorting_id[i]+1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]+1]]+=(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+2]]][i]
                -combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-1]]][i])/
                (combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
            }
            if(init_dis[sorting_id[i][min_sorting_id[i]-1]]!=numeric_limits<double>::infinity()){
                init_dis[sorting_id[i][min_sorting_id[i]-1]]+=(combine_fit[this_rank[sorting_id[i][min_sorting_id[i]+1]]][i]
                -combine_fit[this_rank[sorting_id[i][min_sorting_id[i]-2]]][i])/
                (combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
            }
        }
       
        front_size--;
        for(int i=0;i<obj_num;i++)
            sorting_id[i].erase(sorting_id[i].begin()+min_sorting_id[i]);
        for(int i=0;i<obj_num;i++)
            for(int j=0;j<front_size;j++)
                if(sorting_id[i][j]>min_sol_pos)
                    sorting_id[i][j]--;

        this_rank.erase(this_rank.begin()+min_sol_pos);
        init_dis.erase(init_dis.begin()+min_sol_pos);
        
    }
    
    archive=dd2(rank_0.size());
    d1 tmp_region_id=combine_region_id;
    // archive_region_id=d1(rank_0.size());
    combine_region_id=d1(rank_0.size());
    for(int i=0;i<rank_0.size();i++){
        archive[i]=combine_sol[rank_0[i]];
        combine_region_id[i]=tmp_region_id[rank_0[i]];
        // archive_region_id[i]=tmp_region_id[rank_0[i]];
    }
    
    int deleted_sol_size=deleted_sol.size();
    dd1 deleted_init_dis(deleted_sol_size,0.0);
    for(int i=0;i<deleted_sol_size;i++)
        for(int k=0;k<obj_num;k++)
            for(int j=0;j<front_size;j++)
                if(combine_fit[deleted_sol[i]][k]<combine_fit[this_rank[sorting_id[k][j]]][k]){
                    deleted_init_dis[i]+=(combine_fit[this_rank[sorting_id[k][j]]][k]
                                        -combine_fit[this_rank[sorting_id[k][j-1]]][k])/
                                        (combine_fit[this_rank[sorting_id[k][front_size-1]]][k]
                                        -combine_fit[this_rank[sorting_id[k][0]]][k]);
                    break;
                }
        
    quick_sort(deleted_init_dis,deleted_sol,0,deleted_sol_size-1);
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(deleted_sol.begin(),deleted_sol.end());
    reverse(this_rank.begin(),this_rank.end());

    this_rank.insert(this_rank.end(),deleted_sol.begin(),deleted_sol.end());
    
}
double sede::Euclidean_dis(dd1& sol1,dd1& sol2){
    double sum=0.0;
    for(int i=0;i<obj_num;i++)
        sum+=pow(sol1[i]-sol2[i],2.0);
    return sqrt(sum);
}



void sede::run(){
    int combine_time=numIter/region_num;
    int time_count=combine_time;
    double de_pro;
    for(int i=0;i<numRuns;i++){
        
        init();
        
        combine_time=numIter/(region_num-1);
        time_count=combine_time;
        
        arrange();
        
        for(int j=0;j<numIter;j++){
            
            tmp_iter=j;
            
            generate_new_Cr_F();
            
            de_mut();
           
            cal_ev();
            
            select_player();
            
            marketing_survey();
            
            update_uCr_uF();
            
            if(j>time_count && region_num!=1){
                
                combine_region();
                time_count+=combine_time;
               
            }
            
        }
        
    // cout<<archive.size()<<endl;
    // for(int i=0;i<numIter;i++)
    //     cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
    
    // cout<<"archive_size="<<archive.size()<<endl;
    //Debug:max_pareto應改成archive_size
    // for(int i=0;i<max_pareto;i++){
    //     //cout<<"decode"<<i<<endl;
    //     real_sol[i]=decode(archive[i]);
    // }    
    // cout<<archive.size()<<endl;
    dd2 output_archive(max_pareto);
    output_archive=fitness(archive);
    //輸出該RUN的所有解
        ofstream output_obj;
        output_obj.open("pareto/"+func_name+"/sede/"+func_name+"_sede_"+to_string(i)+".txt");
          for(int k=0;k<max_pareto;k++){
            for(int l=0;l<obj_num;l++)
                output_obj<<output_archive[k][l]<<" ";
            // fixed<<setprecision(20)<<
            output_obj<<endl;
        }
        output_obj.close();
    }
    
    
}

void sede::marketing_survey(){
    for (int i = 0; i < region_num; i++)
        if (tb[i] > 1)
            tb[i] = 1.0;
    // combine_region_id=archive_region_id;
    archive_delete();
    // decom_archive_delete();
    
    // int add_sv=0;

    // dd2 tmp_sampleV(searcher_num*region_num*sample_num);
    // for(int i=0;i<searcher_num;i++)
    //     for(int j=0;j<region_num;j++)
    //         for(int k=0;k<sample_num;k++){
    //             tmp_sampleV[add_sv]=sampleV[i][j][k];
    //             add_sv++;
    //         }
    
    // archive_update(tmp_sampleV);
    
    // Euclidean_archive_d(rank_0);//暫定使用Euclidean
}

void sede::quick_sort(dd1& value,d1& this_rank,int left,int right){
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
void sede::quick_sort_obj(int obj_level,d1& obj_id,d1& sorting_id,int left,int right){
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

        pivot = combine_fit[obj_id[right]][obj_level];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (combine_fit[obj_id[j]][obj_level] <= pivot)
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

sede::dd2 sede::fitness(dd2 & all_sol) {
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
       