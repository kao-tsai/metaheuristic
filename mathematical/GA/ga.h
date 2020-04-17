#ifndef __GA_H_INCLUDED__
#define __GA_H_INCLUDED__
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include "../test_problem.h"

using namespace std;
# define EPS 1.0e-14

class ga{
	public:
		typedef vector<int> solution;
        typedef vector<solution> population;
		ga(int,int,int,double,double,string);
		void run();
	private:
		void init(vector<population>&);
		void crossover(vector<population>&);
        void mutation(vector<population>&);
        vector<vector<double>> fitness(vector<vector<double>>&);
        void TS(vector<population>&);
        vector<population> fast_non_dominated(vector<population>&,vector<vector<double>>);
        solution crowding_dis_assign(solution);
        void obj_sorting(vector<double>&,vector<int>&);
        vector<vector<double>> decode(vector<population>);
        void quick_sort(int,solution&,solution&,int,int);
        void quick_sort_dis(vector<double>&,solution&,int,int);

private:
    
	int numRuns;
	int numIter;
	vector<vector<double>> obj_val;

    int population_num;
	double mutation_pro;
    double crossover_pro;
    string func_name;
    int dimension; //每一個解的維度
    int obj_num; //目標數量
    vector<double> upperbound;
    vector<double> lowerbound;
    solution nbits;
    test_problem func;
    std::random_device rd;
    std::default_random_engine gen = std::default_random_engine(rd());
    
};

ga::ga(int xNumRuns,
int xNumIter,
int xPopulationNum,
double xcrossover_pro,
double xmutation_pro,
string xfunc_name
)
{
    
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
    population_num=xPopulationNum;
    func_name=xfunc_name;
    func=test_problem(func_name);
    dimension=func.idimension;
    obj_num=func.iobj_num;
    upperbound=func.ub;
    lowerbound=func.lb;
    nbits=func.xbits;
    crossover_pro=xcrossover_pro;
    // mutation_pro=xmutation_pro;
    mutation_pro=1.0/(nbits[0]*dimension);
    // if(func_name=="SCH")
    //     mutation_pro=0.08;

}
vector<vector<double>> ga::decode(vector<population> all_sol){
    vector<vector<double>> real_sol(all_sol.size(),vector<double>(dimension));
    int j, k;
    int sum;
    for(int i=0;i<all_sol.size();i++)
    {
        for (j=0; j<dimension; j++)
        {
            sum=0;
            for (k=0; k<nbits[j]; k++)
            {
                if (all_sol[i][j][k]==1)
                {
                    sum += pow(2,nbits[j]-1-k);
                }
            }
           real_sol[i][j] = lowerbound[j] + (double)sum*(upperbound[j] - lowerbound[j])/(double)(pow(2,nbits[j])-1);
        }
    }
    return real_sol;
}

void ga::init(vector<population>& currentSol){
    double rnd;
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dimension; j++)  
            for (int k=0; k<nbits[j]; k++)
            {
                // rnd=(double)rand()/RAND_MAX;
                std::uniform_real_distribution<double> dis(0,1);
                rnd=dis(gen);
                if (rnd <= 0.5)
                    currentSol[i][j][k] = 0;              
                else               
                    currentSol[i][j][k] = 1;        
            }
}

void ga::TS(vector<population> &all_sol){
    vector<population> new_pop(population_num);

    int parent1,parent2;
    //tournament size=2
    for(int i=0;i<population_num;i++){

        parent1=rand()%population_num;
        do{
            parent2=rand()%population_num;
        }while(parent1==parent2);

        if(parent1<parent2)
            new_pop[i]=all_sol[parent1];
        else
            new_pop[i]=all_sol[parent2];
    }
    all_sol=new_pop;
}


void ga::crossover(vector<population> &parent){
    //two-point crossover
    /*
    vector<population>child;
    child=parent;
    double rnd;
    int temp, site1, site2;
    int y=0;
  
    for (int i=0; i<population_num/2; i++)
    {
        for(int j=0;j<dimension;j++){
            rnd = (double)rand()/RAND_MAX;
            if (rnd <= crossover_pro)
            {
                site1 = rand()%(nbits[j]-1);
                site2=rand()%(nbits[j]-1);
                if(site1>site2){
                    temp=site1;
                    site1=site2;
                    site2=temp;
                }
                for(int k=site1;k<=site2;k++){
                    child[y][j][k]=parent[y+1][j][k];
                    child[y+1][j][k]=parent[y][j][k];
                }
    
            }
        }
        y=y+2;
        
    }
    parent=child;
    return;*/
    

    //-------------single point crossover----------------//
    vector<population>child;
    child=parent;
    double rnd;
    int temp, site1;
    int y=0;
  
    for (int i=0; i<population_num/2; i++)
    {
        for(int j=0;j<dimension;j++){
            // rnd = (double)rand()/RAND_MAX;
            std::uniform_real_distribution<double> dis(0,1);
            rnd=dis(gen);
            if (rnd <= crossover_pro)
            {
                std::uniform_int_distribution<int> gene(0,nbits[j]-2);
                site1=gene(gen);
                // site1 = rand()%(nbits[j]-1);
                
                for(int k=0;k<=site1;k++){
                    child[y][j][k]=parent[y+1][j][k];
                    child[y+1][j][k]=parent[y][j][k];
                }
    
            }
        }
        y=y+2;
        
    }
    parent=child;
    return;
    
    //不分維度各維度一起crossover
    // vector<population>child;
    // child=parent;
    // double rnd;
    // int temp, site1;
    // int y=0;
    // int chrom=0,c,w_dim,w_bit;
    // for(int i=0;i<dimension;i++)
    //     chrom+=nbits[i];
    // for (int i=0; i<population_num/2; i++)
    // {

    //     rnd = (double)rand()/RAND_MAX;
    //     if (rnd <= crossover_pro)
    //     {
    //         rnd=(double)rand()/RAND_MAX;
    //         c=floor(rnd*(chrom+10));
            
    //         site1 = c;
    //         if(site1>=chrom)
    //             site1=site1/2;
            
    //         for(int j=0;j< chrom;j++){
                
    //             w_dim=j/nbits[0];           
    //             w_bit=j%nbits[0];
    //             if(j>site1-1){    
    //                 child[y][w_dim][w_bit]=parent[y+1][w_dim][w_bit];        
    //                 child[y+1][w_dim][w_bit]=parent[y][w_dim][w_bit];
    //             }
    //         }
    //     }
    //     y=y+2;  
    // }
    // parent=child;
    // return;
    

    // vector<population>child;
    // child=parent;
    // double rnd;
    // int temp, site1;
    // int y=0;
  
    // for (int i=0; i<population_num/2; i++)
    // {
    //     std::uniform_real_distribution<double> dis(0,1);
    //     rnd=dis(gen);
    //     // rnd = (double)rand()/RAND_MAX;
    //     if (rnd <= crossover_pro)
    //     {
    //         for(int j=0;j<dimension;j++){

    //         // site1 = rand()%(nbits[j]-1);
    //         std::uniform_int_distribution<int> gene(0,nbits[j]-2);
    //         site1=gene(gen);
    //             for(int k=0;k<=site1;k++){
    //                 child[y][j][k]=parent[y+1][j][k];
    //                 child[y+1][j][k]=parent[y][j][k];
    //             }
    //         }
    //     }
        
    //     y=y+2;
        
    // }
    // parent=child;
    // return;
}
void ga::mutation(vector<population>& all_sol){

    double prob;
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dimension; j++)
            for (int k=0; k<nbits[j]; k++)
            {
                std::uniform_real_distribution<double> dis(0,1);
                prob=dis(gen);
                // prob = (double)rand()/RAND_MAX;
                if (prob <=mutation_pro)
                    all_sol[i][j][k]=(all_sol[i][j][k]+1)%2;            
            }

    return;

}
//-----------------------------------------------------run nsga2 algorithm----------------------------------------------------------//
void ga::run(){
    
    vector<double>iter_obj_avg(numIter,0.0);
    vector<population> old_binary_sol(population_num,population(dimension));
    vector<population> new_binary_sol;
    vector<vector<double>> real_sol(population_num,vector<double>(dimension));
    
   for(int i=0;i<population_num;i++)
       for(int j=0;j<dimension;j++)     
           old_binary_sol[i][j].resize(nbits[j]);
    
    for(int i=0;i<numRuns;i++){
        
        init(old_binary_sol);
        
        real_sol=decode(old_binary_sol);
        //-----------------------------------------------//
        //印出所有人口測試
/*
        for(int a=0;a<population_num;a++){
            for(int b=0;b<dimension;b++)
                cout<<real_sol[a][b]<<" ";
            cout<<endl;
        }*/
         //-----------------------------------------------//
        new_binary_sol=fast_non_dominated(old_binary_sol,real_sol);
        
        for(int j=0;j<numIter;j++){
            
            TS(new_binary_sol);
            
            crossover(new_binary_sol);
            
            mutation(new_binary_sol);
    
            //IGD
            new_binary_sol.insert(new_binary_sol.end(),old_binary_sol.begin(),old_binary_sol.end());//combination
            real_sol=decode(new_binary_sol);
            old_binary_sol=fast_non_dominated(new_binary_sol,real_sol);
            
            new_binary_sol=old_binary_sol;
            //iter_obj_avg[j]+=best_obj_value;
        }
        real_sol=decode(new_binary_sol);
        obj_val=fitness(real_sol);
        //輸出該RUN的所有解
        ofstream output_obj;
        output_obj.open("pareto/"+func_name+"/bnsga2/"+func_name+"_bnsga2_"+to_string(i)+".txt");
          for(int k=0;k<population_num;k++){
            for(int l=0;l<obj_num;l++)
                output_obj<<obj_val[k][l]<<" ";
            output_obj<<endl;
        }
        output_obj.close();
    }
    /*
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
    */
    // return obj_val;

}
//------------------------------------------------------------------------------------------------------------------------------------------//
vector<ga::population> ga::fast_non_dominated(vector<population>& new_b_sol,vector<vector<double>> combine_sol){
    /*
    combine_sol[0][0]=1;
    combine_sol[1][0]=-1;
    combine_sol[2][0]=2;
    combine_sol[3][0]=0.5;
    combine_sol[4][0]=3;*/
    //帶入公式算出每個目標的值
    obj_val=fitness(combine_sol);
    int combine_num=combine_sol.size();
    
    //印出所有目標的值(記住要刪除)//
    /*
    for(int i=0;i<combine_num;i++){
        for(int j=0;j<obj_num;j++)
            cout<<obj_val[i][j]<<" ";
        cout<<endl;
    }
    */
    //---------------------------------------//

    //紀錄dominated數量和dominate set
    vector<vector<int>> set(combine_num,vector<int>(1,0));
    for(int i=0;i<combine_num-1;i++)
    {
        for(int j=i+1;j<combine_num;j++)
        
            if(obj_val[i][0]<=obj_val[j][0]&&obj_val[i][1]<=obj_val[j][1]){
                if(obj_val[i][0]==obj_val[j][0]&&obj_val[i][1]==obj_val[j][1])
                    continue;
                else{
                    set[i].push_back(j);
                    set[j][0]++;
                }
            }
            else if(obj_val[i][0]>=obj_val[j][0]&&obj_val[i][1]>=obj_val[j][1]){
                set[i][0]++;
                set[j].push_back(i);
            }
    }
   
    //(記住要山)----------------------//
/*
    for(int i=0;i<set.size();i++){
        for(int j=0;j<set[i].size();j++)
            cout<<set[i][j]<<" ";
        cout<<endl;
    }
*/
    //-----------------------------------//
    int new_population=0;
    //開始排rank
    population rank;
    solution temp_rank;
    while(new_population<population_num){
        for(int i=0;i<combine_num;i++)       
            if(set[i][0]==0){
                temp_rank.push_back(i);
                set[i][0]=-1;
            }
        rank.push_back(temp_rank);
        for(int i=0;i<temp_rank.size();i++)
            for(int j=1;j<set[temp_rank[i]].size();j++)
                set[set[temp_rank[i]][j]][0]--;
        new_population+=temp_rank.size();
        temp_rank.clear();
    }
    
    //(記住要刪除)--------------------//
    /*
    cout<<"here is rank:"<<endl;
    for(int i=0;i<rank.size();i++){
        for(int j=0;j<rank[i].size();j++)
            cout<<rank[i][j]<<" ";
        cout<<endl;
    }*/
    
    //-----------------------------------//
    /*
    cout<<rank.size()<<endl;
    if(rank.size()==6)
    for(int i=0;i<6;i++)
        cout<<obj_val[rank[i][0]][0]<<" "<<obj_val[rank[i][0]][1]<<endl;
    */
    //做crowding_distance
    for(int i=0;i<rank.size();i++){
        if(rank[i].size()>1)
            rank[i]=crowding_dis_assign(rank[i]);
    }
    
    //(記住要刪除)--------------------//
    /*
    cout<<"here is rank after crowding distance:"<<endl;
    for(int i=0;i<rank.size();i++){
        for(int j=0;j<rank[i].size();j++)
            cout<<rank[i][j]<<" ";
        cout<<endl;
    }*/
    //-----------------------------------//
    //清除多餘的人口
    while(new_population>population_num){
        //cout<<"清除中"<<endl;
        rank[rank.size()-1].pop_back();
        new_population--;
    }
     
    //最終排序人口回傳
    vector<population> final_population(population_num);
    int x=0;
    for(int i=0;i<rank.size();i++){
        for(int j=0;j<rank[i].size();j++){
            final_population[x]=new_b_sol[rank[i][j]];
            x++;
        }
    }
    
    //記住要刪除---------------------------------//
    /*for(int i=0;i<final_population.size();i++)
        for(int j=0;j<dimension;j++)
            cout<<final_population[i][j]<<endl;*/
    //-----------------------------------------------//
    return final_population;
}



vector<vector<double>> bnsgaii::fitness(vector<vector<double>> & all_sol) {
    
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

//--------------------------------------------Use Quick Sort------------------------------------------------------//

void bnsgaii::quick_sort(int obj_level,solution& obj_id,solution& sorting_id,int left,int right){
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

        pivot = obj_val[obj_id[right]][obj_level];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (obj_val[obj_id[j]][obj_level] <= pivot)
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

        quick_sort(obj_level,obj_id,sorting_id, left, index-1);
        quick_sort(obj_level,obj_id,sorting_id, index+1, right);
    }
}


void bnsgaii::quick_sort_dis(vector<double>& init_dis,solution& this_rank,int left,int right){
    int index;
    int int_tmp;
    double double_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        double_tmp = init_dis[right];
        init_dis[right] = init_dis[index];
        init_dis[index] = double_tmp;

        int_tmp = this_rank[right];
        this_rank[right] = this_rank[index];
        this_rank[index] = int_tmp;

        pivot = init_dis[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (init_dis[j] <= pivot)
            {
                i+=1;
                double_tmp = init_dis[j];
                init_dis[j] = init_dis[i];
                init_dis[i] = double_tmp;
                
                int_tmp = this_rank[j];
                this_rank[j] = this_rank[i];
                this_rank[i] = int_tmp;
            }
        }
        index=i+1;
        double_tmp = init_dis[index];
        init_dis[index] = init_dis[right];
        init_dis[right] = double_tmp;

        int_tmp = this_rank[index];
        this_rank[index] = this_rank[right];
        this_rank[right] = int_tmp;

        quick_sort_dis(init_dis,this_rank, left, index-1);
        quick_sort_dis(init_dis,this_rank, index+1, right);
    }
}

#endif
       