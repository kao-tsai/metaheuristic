#ifndef __SE_H_INCLUDED__
#define __SE_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include<iomanip>
#include"../test_problem.h"
using namespace std;

class se{
	public:
		typedef vector<int> d1;
        typedef vector<d1> d2;
        typedef vector<d2> d3;
        typedef vector<d3> d4;
        typedef vector<double> dd1;
        typedef vector<dd1> dd2;
        typedef vector<dd2> dd3;
        typedef vector<dd3> dd4;
		se(int,int,int,int,int,int,int,string,double,double);
		void run();
	private:
		void init();
		
        void ga_crossover();
        void ga_mut();
        void arrange();
        void cal_ev();
        void choose_sample(int);
        void select_player();
        void marketing_survey();
        void init_region_bit();
        dd1 decode(d2&);
        void cal_cdr_fit();
        dd2 fitness(dd2&);
        void crowding_dis_assign(d1&);
        void quick_sort(dd1&,d1&,int,int);
        void quick_sort_obj(int,d1&,d1&,int,int);
        void archive_delete();
private:
	int numRuns;
	int numIter;
	int dim;

    d3 searcher;
    d4 sample;
    vector<d4> sampleV;
    d2 region_bit;
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

    d3 archive;//Best Non-Dominated Solutions

    d3 combine_sol;//暫存Sample+SampleV+Archive
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

    string func_name;
    test_problem func;
};

se::se(int xNumRuns,
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
    
    region_num=xRegion;
    searcher_num=in_searcher_num;
    sample_num=in_sample_num;
    tour_size=in_tournament_size;
    func_name=in_func_name;
    func=test_problem(func_name);
    // func_init(func_name);
    dim=func.idimension;
    obj_num=func.iobj_num;
    max_pareto=in_max_pareto;
    nbits=func.xbits;
    lowerbound=func.lb;
    upperbound=func.ub;
    full_sample_num=sample_num*region_num;
    tatb_pressure=2;
    dominate_pressure=3;
    crossover_pro=in_crossover_pro;
    // mutation_pro=in_mutation_pro;
    mutation_pro=1.0/(nbits[0]*dim);
}


void se::init(){

    tb.assign(region_num,0.0);
    ta.assign(region_num,1.0);
    searcher_belong.assign(searcher_num,0);
    
    searcher.assign(searcher_num,d2(dim,d1(nbits[0],0)));
    
    sample.assign(region_num,d3(sample_num,d2(dim,d1(nbits[0],0))));
    sampleV.assign(searcher_num, d4(region_num, d3(sample_num*2, d2(dim, d1(nbits[0],0)))));
    sample_cdr_fit.assign(region_num,dd1(sample_num,0.0));
    sampleV_cdr_fit.assign(searcher_num,dd2(region_num,dd1(sample_num*2,0.0)));
    region_best_fit.assign(region_num,0);

    
    archive=d3(searcher_num+region_num*sample_num);
	for (int i = 0; i < searcher_num; i++){
        for(int j=0;j<dim;j++)
            for(int k=0;k<nbits[j];k++)
                searcher[i][j][k]=rand()%2;
        archive[i]=searcher[i];
    }
    
}
se::dd1 se::decode(d2& sol){
    dd1 real_sol(dim,0.0);
    int j, k;
    int sum;

    for (j=0; j<dim; j++)
    {
        sum=0;
        for (k=0; k<nbits[j]; k++)
            if (sol[j][k]==1)
                sum += pow(2,nbits[j]-1-k);
    
        real_sol[j] = lowerbound[j] + (double)sum*(upperbound[j] - lowerbound[j])/(double)(pow(2.0,nbits[j])-1);
        
    }

    return real_sol;
}
void se::init_region_bit(){
    region_bit.assign(region_num,d1(clip_bit_num,0));

    for(int i=0;i<region_num;i++){
        int n=clip_bit_num-1;
        int tmp=i;
        while(n>=0){
            region_bit[i][n]=tmp%2;
            tmp/=2;
            n--;
        }
    }
    //記得刪除----------------//
    // for(int i=0;i<region_num;i++){
    //     for(int j=0;j<clip_bit_num;j++)
    //         cout<<region_bit[i][j]<<endl;
    //     cout<<endl;
    // }

}
void se::arrange(){
    clip_bit_num=log2(region_num);
    
    init_region_bit();
    //--------------------專注於第1個維度切割---------------------//
    for(int i=0;i<searcher_num;i++){
        searcher_belong[i]=i%region_num;
        for(int j=0;j<clip_bit_num;j++)
            searcher[i][0][j]=region_bit[searcher_belong[i]][j];
        //cout<<"Searcher "<<i<<" : "<<searcher_belong[i]<<endl;
    }
    
    //-------------------平均切割每個維度-------------------------//
    /*
    for(int i=0;i<searcher_num;i++){
        searcher_belong[i]=i%region_num;
        for(int j=0;j<dim;j++)
            searcher[i][j][0]=region_bit[searcher_belong[i]][j];
    }*/
    //----------------------------------------------------------//
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            for(int k=0;k<dim;k++)
                for(int l=0;l<nbits[k];l++)
                    if(l<clip_bit_num && k==0)
                        sample[i][j][k][l]=region_bit[i][l];
                    else
                        sample[i][j][k][l]=rand()%2;
            archive[searcher_num+i*sample_num+j]=sample[i][j];
        }
    for(int i=0;i<searcher_num;i++){
        tb[searcher_belong[i]]++;
        ta[searcher_belong[i]]=1.0;
    }
                
    EV.assign(searcher_num,dd1(region_num,0.0));

   //取得所有二元解的實數解並計算Objective Value
    int archive_size=archive.size();
    dd2 archive_real_sol(archive_size);
    archive_fit.clear();//可能不需要
    for(int i=0;i<archive_size;i++)
        archive_real_sol[i]=decode(archive[i]);
    archive_fit=fitness(archive_real_sol);
    

    //紀錄dominated數量和dominate set
    d2 set(archive_size,d1(1,0));
    for(int i=0;i<archive_size-1;i++) 
        for(int j=i+1;j<archive_size;j++)
            if(archive_fit[i][0]<=archive_fit[j][0]&&archive_fit[i][1]<=archive_fit[j][1])
                if(archive_fit[i][0]==archive_fit[j][0]&&archive_fit[i][1]==archive_fit[j][1])
                    continue;
                else{
                    set[i].push_back(j);
                    set[j][0]++;
                }
            else if(archive_fit[i][0]>=archive_fit[j][0]&&archive_fit[i][1]>=archive_fit[j][1]){
                set[i][0]++;
                set[j].push_back(i);
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
    rank_0=rank[0];
    d3 tmp_archive;
    tmp_archive=archive;
    archive=d3(rank_0.size());
     for(int i=0;i<rank_0.size();i++)
            archive[i]=tmp_archive[rank_0[i]];
    
}


void se::ga_crossover(){
    double rnd;
    int crossover_point,s_pos;
    for(int i=0;i<searcher_num;i++)
        for(int j=0;j<region_num;j++)
            for(int k=0;k<sample_num;k++){
                for(int l=0;l<dim;l++){
                    rnd=(double)rand()/RAND_MAX;
                    if(rnd<crossover_pro){
                        crossover_point=rand()%(nbits[l]-1)+1;
                        //若維度在第1維度要固定切割的區域
                        if(l==0){
                            for (int m = 0; m < clip_bit_num; m++) {
                                sampleV[i][j][k*2][l][m] = region_bit[j][m];
                                sampleV[i][j][k*2+1][l][m] = region_bit[j][m];
                            }
                            for (int m = clip_bit_num; m < nbits[l]; m++) 
                                if (l < crossover_point) {
                                    sampleV[i][j][k*2][l][m] = searcher[i][l][m];
                                    sampleV[i][j][k*2+1][l][m] = sample[j][k][l][m];
                                }
                                else {
                                    sampleV[i][j][k*2][l][m] = sample[j][k][l][m];
                                    sampleV[i][j][k*2+1][l][m] = searcher[i][l][m];
                                }
                        }
                        //其他維度正常crossover
                        else{
                            for (int m = 0; m < nbits[l]; m++) {
                                if (l < crossover_point) {
                                    sampleV[i][j][k*2][l][m] = searcher[i][l][m];
                                    sampleV[i][j][k*2+1][l][m] = sample[j][k][l][m];
                                }
                                else {
                                    sampleV[i][j][k*2][l][m] = sample[j][k][l][m];
                                    sampleV[i][j][k*2+1][l][m] = searcher[i][l][m];
                                }
                            }
                        }
                    }
                    //rnd>crossover probability 不做Crossover
                    else{
                        sampleV[i][j][k*2]=searcher[i];
                        sampleV[i][j][k*2+1]=sample[j][k];
                        for(int m=0;m<clip_bit_num;m++)
                            sampleV[i][j][k*2][0][m]=region_bit[j][m];
                        break;
                    }
                }
            }
}
void se::ga_mut(){
    double rnd;
    for (int i = 0; i < searcher_num; i++) 
        for (int j = 0; j < region_num; j++) 
            for (int k = 0; k < 2*sample_num; k++) 
                for(int l=0;l<dim;l++)
                    for(int m=0;m<nbits[l];m++){
                        rnd=(double)rand()/RAND_MAX;
                        if(rnd<mutation_pro){
                            if(l==0 && m<clip_bit_num)
                                continue; 
                            sampleV[i][j][k][l][m]=!sampleV[i][j][k][l][m];
                        }
                    }
}

//考慮三種計算CDR的方法
//1.將所有sample和sampleV加入archive並計算crowding_distance
/*2.每個點依次放入archive判斷單點在archive內的優異度，
以加入後的新擁擠距離與最小擁擠距離相減並加上新點之分佈平均*/
//3.利用Spread()公式判斷各點加入archive後的優異度
void se::cal_cdr_fit(){

    combine_sol=archive;
    int archive_size=archive.size();
    d1 store_sample_pos;
    int sample_end=archive_size-1;
    
    //暫時將Sample加入archive內比較其優異度，若已有相同解在內則無須加入，紀錄位置即可
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            bool flag=true;
            for(int k=0;k<archive_size;k++)
                if(combine_sol[k]==sample[i][j]){
                    store_sample_pos.push_back(k);
                    flag=false;
                    break;
                } 
            if(flag){
                combine_sol.push_back(sample[i][j]);
                sample_end++;
                store_sample_pos.push_back(sample_end);
            }
        }
    
    //紀錄所有combine_sol的region_id
    combine_region_id=d1(combine_sol.size(),-1);
    for(int i=0;i<store_sample_pos.size();i++){
        combine_region_id[store_sample_pos[i]]=i/sample_num;
        // cout<<combine_region_id[store_sample_pos[i]]<<endl;
    }
    for(int i=0;i<combine_sol.size();i++){
        if(combine_region_id[i]==-1)
            for(int j=0;j<region_num;j++){
                // cout<<combine_sol[i][0][0]<<combine_sol[i][0][1]<<endl;
                if(equal(combine_sol[i][0].begin(),combine_sol[i][0].begin()+clip_bit_num,region_bit[j].begin())){
                    combine_region_id[i]=j;
                    break;
                }
            }
        
    }
    
    // cout<<"combine_sol_size1:"<<combine_sol.size()<<endl;
    //將sampleV暫時加入archive內比較他的優異度
    int sampleV_start=sample_end+1;
    // cout<<"SampleV_start:"<<sampleV_start<<endl;
    for(int i=0;i<region_num;i++)
        for(int j=0;j<searcher_num;j++)
            for(int k=0;k<sample_num;k++){
                combine_sol.push_back(sampleV[j][i][2*k]);
                // combine_region_id[sampleV_start]=i;
                combine_region_id.push_back(i);
                combine_sol.push_back(sampleV[j][i][2*k+1]);
                // combine_region_id[sampleV_start+1]=i;
                combine_region_id.push_back(i);
              
            }
    sampleV_start=sample_end+1;
    // cout<<"combine_sol_size2:"<<combine_sol.size()<<endl;
    

    //取的所有二元解的實數解並取得Objective Value
    int combine_size=combine_sol.size();
    dd2 combine_real_sol(combine_size);
    combine_fit.clear();//可能不需要
    for(int i=0;i<combine_size;i++)
        combine_real_sol[i]=decode(combine_sol[i]);
    combine_fit=fitness(combine_real_sol);
    // cout<<"combine_fit_size:"<<combine_fit.size()<<endl;

    //紀錄dominated數量和dominate set
    d2 set(combine_size,d1(1,0));
    for(int i=0;i<combine_size-1;i++) 
        for(int j=i+1;j<combine_size;j++)
            if(combine_fit[i][0]<=combine_fit[j][0]&&combine_fit[i][1]<=combine_fit[j][1])
                if(combine_fit[i][0]==combine_fit[j][0]&&combine_fit[i][1]==combine_fit[j][1])
                    continue;
                else{
                    set[i].push_back(j);
                    set[j][0]++;
                }
            else if(combine_fit[i][0]>=combine_fit[j][0]&&combine_fit[i][1]>=combine_fit[j][1]){
                set[i][0]++;
                set[j].push_back(i);
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

    //做crowding_distance
    for(int i=0;i<rank.size();i++)
        if(rank[i].size()>2)
            crowding_dis_assign(rank[i]);
    
    //1.嘗試填滿archive(bugfix:2020/03/18)
    if(rank[0].size()>max_pareto)
        rank_0=rank[0];
    else{
        rank_0=rank[0];
        int insert_rank_i=1;
        while(rank_0.size()<max_pareto && insert_rank_i<rank.size()){
            rank_0.insert(rank_0.end(),rank[insert_rank_i].begin(),rank[insert_rank_i].end());
            insert_rank_i++;
        }
        while(rank_0.size()>max_pareto)
            rank_0.pop_back();
    }
    //2.只存rank1的解但可能會出錯
    //rank_0=rank[0];

    ///----------------------------------------------------///
    // cout<<"rank size "<<rank.size()<<endl;
    combine_region_id_2d=d2(region_num);
    for(int i=0;i<rank.size();i++)
        for(int j=0;j<rank[i].size();j++)
            combine_region_id_2d[combine_region_id[rank[i][j]]].push_back(rank[i][j]);
        
   
    //計算所有解的CDR
    double giving_level=(double)combine_size;
    double rank_size=(double)rank.size();
    combine_cdr_fit.assign(combine_size,0.0);
    for(int i=0;i<rank_size;i++)
        for(int j=0;j<rank[i].size();j++){
            combine_cdr_fit[rank[i][j]]=giving_level*pow(dominate_pressure,(rank_size-i));
            giving_level=giving_level-1;
        }
    
    //存入Sample_cdr_fit和SampleV_cdr_fit
    int count=0;
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++){
            sample_cdr_fit[i][j]=combine_cdr_fit[store_sample_pos[count]];
            count++;
        }
    count=sampleV_start;

    // cout<<"sampleV_start:"<<count<<endl;

    for(int i=0;i<region_num;i++)
        for(int j=0;j<searcher_num;j++)
            for(int k=0;k<sample_num;k++){
                sampleV_cdr_fit[j][i][k*2]=combine_cdr_fit[count];
                sampleV_cdr_fit[j][i][k*2+1]=combine_cdr_fit[count+1];
                count+=2;
            }
    
}

//記得刪除--------------------------------------------------------------------//(試試是否正確)
void se::crowding_dis_assign(d1& this_rank){
   
    int front_size=this_rank.size();
    dd1 init_dis(front_size,0.0);

    d2 sorting_id(obj_num,d1(front_size));
    d2 obj_id(obj_num);
    for(int i=0;i<obj_num;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    }
    
    //計算擁擠程度
    for (int i=0; i<obj_num; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
    
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<front_size-1;j++){
            if(init_dis[j]!=numeric_limits<double>::infinity()){
                if(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]== combine_fit[this_rank[sorting_id[i][0]]][i])
                    init_dis[sorting_id[i][j]]+=0.0;
                else
                    init_dis[sorting_id[i][j]]+=(combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i])/(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
            }
        }
    //距離密度值進行正規化
    for (int j=0; j<front_size; j++){
        if (init_dis[j]!= numeric_limits<double>::infinity())
           init_dis[j] = init_dis[j]/obj_num;
        //    cout<<"init_dis "<<j<<": "<<init_dis[j]<<endl;
    }
    
    quick_sort(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end());
    
}
//可能會刪除
/*
void se::crowding_dis_assign(d1& this_rank){
   
    int solution_num=this_rank.size();
    dd1 init_dis(solution_num,0.0);
    dd2 this_obj_val(solution_num);

    for(int i=0;i<solution_num;i++){
        this_obj_val[i]=combine_fit[this_rank[i]];
        
    }

    d1 tmp_id(solution_num);
    d2 sorting_id(obj_num);
    for(int i=0;i<solution_num;i++)
        tmp_id[i]=i;
    for(int i=0;i<obj_num;i++)
        sorting_id[i]=tmp_id;
    //列行互換
    vector<vector<double>> row_obj_val(obj_num,vector<double>(solution_num));
    for(int i=0;i<obj_num;i++)
        for(int j=0;j<solution_num;j++)
            row_obj_val[i][j]=this_obj_val[j][i];

    //對各目標進行排序
    
    // obj_sorting(row_obj_val[0],sorting_id[0]);
    // sorting_id[1]=sorting_id[0];
    // reverse(sorting_id[1].begin(),sorting_id[1].end());
    // vector<double> origin_obj1;
    // origin_obj1=row_obj_val[1];
    // for(int i=0;i<solution_num;i++)
    //     row_obj_val[1][i]=origin_obj1[sorting_id[1][i]];
    for(int i=0;i<obj_num;i++)
        quick_sort(row_obj_val[i],sorting_id[i],0,solution_num-1);

    //計算擁擠程度
    for (int i=0; i<obj_num; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();

    for(int i=0;i<obj_num;i++)
        for(int j=1;j<solution_num-1;j++){
             if(init_dis[j]!=numeric_limits<double>::infinity()){
                 if(row_obj_val[i][solution_num-1]==row_obj_val[i][0])
                    init_dis[sorting_id[i][j]]+=0.0;
                else
                    init_dis[sorting_id[i][j]]+=(row_obj_val[i][j+1]-row_obj_val[i][j-1])/(row_obj_val[i][solution_num-1]-row_obj_val[i][0]);
             }
        }
    //距離密度值進行正規化
    for (int j=0; j<solution_num; j++){
        if (init_dis[j]!= numeric_limits<double>::infinity())
           init_dis[j] = init_dis[j]/obj_num;
        //記得刪
        cout<<"init_dis "<<j<<" : "<<init_dis[j]<<endl;
    }
    
    quick_sort(init_dis,this_rank,0,solution_num-1);
    reverse(this_rank.begin(),this_rank.end());
}*/
void se::cal_ev(){

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
                Vij[i][j] += sampleV_cdr_fit[i][j][k*2] + sampleV_cdr_fit[i][j][k*2+1];
            }
            Vij[i][j] /= 2*sample_num; 
        }
    }
    //計算所有Searcher各區域的期望值
     for (int i = 0; i < searcher_num; i++) 
        for (int j = 0; j < region_num; j++) 
            EV[i][j] = (tb[j]/ta[j])*Vij[i][j]*(region_best_fit[j]/sample_fit_sum);
        
    // for(int i=0;i<searcher_num;i++)
    //     searcher_fit[i]=evaluation(searcher[i]);
    choose_sample(1);
}
void se::choose_sample(int choose_method){
    //考慮三種作法:
    //1.sample與對應sampleV比較並取最大值
    //2.同區域內取前幾大(sample數)做比較並替換原sample
    //3.採用Tournament Selection來選Sample
    //以下採用(1)
    if(choose_method==1){
        for (int i = 0; i < searcher_num ; i++) 
            for(int j = 0; j < region_num; j++) 
                for (int k = 0; k < sample_num; k++) {
                    if (sampleV_cdr_fit[i][j][k*2] > sample_cdr_fit[j][k]) {
                        // for (int l = clip_bit_num; l < numPatterns; l++)
                        sample[j][k]= sampleV[i][j][k*2];
                        sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k*2];
                    }
                    if (sampleV_cdr_fit[i][j][k*2+1] > sample_cdr_fit[j][k]) {
                        // for (int l = clip_bit_num; l < numPatterns; l++)
                        sample[j][k] = sampleV[i][j][k*2+1];
                        sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k*2+1];
                    }
                }
    }
    //以下採用(2)(未寫)
    else if(choose_method==2){

    }
    //以下採用(3)(2020/03/17 add in)
    else if(choose_method==3){
        for (int i = 0; i < searcher_num ; i++) 
            for(int j = 0; j < region_num; j++) 
                for (int k = 0; k < sample_num; k++) {
                    if (sampleV_cdr_fit[i][j][k*2] > sample_cdr_fit[j][k]) {
                   
                        sample[j][k]= sampleV[i][j][k*2];
                        sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k*2];
                    }
                    if (sampleV_cdr_fit[i][j][k*2+1] > sample_cdr_fit[j][k]) {
                        
                        sample[j][k] = sampleV[i][j][k*2+1];
                        sample_cdr_fit[j][k] = sampleV_cdr_fit[i][j][k*2+1];
                    }
                }
    }
    
}
void se::select_player(){
    for (int i = 0; i < region_num; i++)
        tb[i]++;

    for (int i = 0; i < searcher_num; i++) {

        int first_select = rand() % region_num;
        double first_ev = EV[i][first_select];
        
      
        int second_select= rand() % region_num;
        if (EV[i][second_select] >  first_select) {
            first_select = second_select;
            first_ev = EV[i][first_select];
        }
       
        //決定每個Searcher，有兩種作法
        //1.選取最好------->(1)選前幾好的(X)，(2)一律選最好的(V)
        //2.使用Tournament Selection選取
        //以下採用(1.1)
        if(combine_region_id_2d[first_select].size()>0){
            
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
            searcher[i]=sample[first_select][best_pos];
        }
         
        // update region_it[i] and region_hl[i];
        ta[first_select]++;
        tb[first_select] = 1;
    }
}
void se::archive_delete(){

       while(rank_0.size()>max_pareto){
            crowding_dis_assign(rank_0);
            rank_0.pop_back();

        }
        archive=d3(rank_0.size());
        for(int i=0;i<rank_0.size();i++)
            archive[i]=combine_sol[rank_0[i]];

}
void se::run(){
    
    dd1 iter_obj_avg(numIter,0.0);
    for(int i=0;i<numRuns;i++){
        
        init();
        
        arrange();
       
        for(int j=0;j<numIter;j++){
           
            ga_crossover();
            
            ga_mut();
            // cout<<j<<" : hello"<<endl;
            cal_ev();
            
            select_player();
            
            marketing_survey();
            
            // iter_obj_avg[j]+=best_obj_val;
        }

    // for(int i=0;i<numIter;i++)
    //     cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
    dd2 real_sol(max_pareto);
    // cout<<"archive_size="<<archive.size()<<endl;
    //Debug:max_pareto應改成archive_size
    for(int i=0;i<max_pareto;i++){
        // cout<<"decode"<<i<<endl;
        real_sol[i]=decode(archive[i]);
    }    

    // cout<<"decode success"<<endl;

    dd2 output_archive(max_pareto);
    output_archive=fitness(real_sol);
    //輸出該RUN的所有解
        ofstream output_obj;
        output_obj.open("pareto/"+func_name+"/se/"+func_name+"_se_"+to_string(i)+".txt");
          for(int k=0;k<max_pareto;k++){
            for(int l=0;l<obj_num;l++)
                output_obj<<output_archive[k][l]<<" ";
            output_obj<<endl;
        }
        output_obj.close();
    }
}

void se::marketing_survey(){
    for (int i = 0; i < region_num; i++)
        if (tb[i] > 1)
            tb[i] = 1.0;

    archive_delete();
}

void se::quick_sort(dd1& value,d1& this_rank,int left,int right){
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
void se::quick_sort_obj(int obj_level,d1& obj_id,d1& sorting_id,int left,int right){
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

se::dd2 se::fitness(dd2 & all_sol) {
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
}
#endif
       