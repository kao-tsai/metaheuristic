#ifndef __BNSGAII_H_INCLUDED__
#define __BNSGAII_H_INCLUDED__
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>

#include "dataset.cpp"
#include"function_init.cpp"


using namespace std;
# define EPS 1.0e-14

class bnsgaii{
	public:
		typedef vector<int> solution;
        typedef vector<solution> population;
		bnsgaii(int,int,int,double,double,string);
		vector<vector<double>> run();
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
    string problem_func;
    int dimension; //每一個解的維度
    int obj_num; //目標數量
    vector<double> upperbound;
    vector<double> lowerbound;
    solution nbits;
    
};

bnsgaii::bnsgaii(int xNumRuns,
int xNumIter,
int xPopulationNum,
double xcrossover_pro,
double xmutation_pro,
string xproblem_func
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
    crossover_pro=xcrossover_pro;
    mutation_pro=xmutation_pro;
    // mutation_pro=1/30;
    population_num=xPopulationNum;
    problem_func=xproblem_func;
    func_init(problem_func);
    dimension=idimension;
    obj_num=iobj_num;
    upperbound=ub;
    lowerbound=lb;
    nbits=xbits;

}
vector<vector<double>> bnsgaii::decode(vector<population> all_sol){
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

void bnsgaii::init(vector<population>& currentSol){
    double rnd;
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dimension; j++)  
            for (int k=0; k<nbits[j]; k++)
            {
                rnd=(double)rand()/RAND_MAX;
                if (rnd <= 0.5)
                    currentSol[i][j][k] = 0;              
                else               
                    currentSol[i][j][k] = 1;        
            }
}

void bnsgaii::TS(vector<population> &all_sol){
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


void bnsgaii::crossover(vector<population> &parent){
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
                for(int k=site1;k<=site1;k++){
                    child[y][j][k]=parent[y+1][j][k];
                    child[y+1][j][k]=parent[y][j][k];
                }
    
            }
        }
        y=y+2;
        
    }
    parent=child;
    return;
}
void bnsgaii::mutation(vector<population>& all_sol){

    double prob;
    for(int i=0;i<population_num;i++)
        for (int j=0; j<dimension; j++)
            for (int k=0; k<nbits[j]; k++)
            {
                prob = (double)rand()/RAND_MAX;
                if (prob <=mutation_pro)
                    all_sol[i][j][k]=(all_sol[i][j][k]+1)%2;            
            }

    return;

}
//-----------------------------------------------------run nsga2 algorithm----------------------------------------------------------//
vector<vector<double>> bnsgaii::run(){
    
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
    }
    /*
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
    */
    return obj_val;

}
//------------------------------------------------------------------------------------------------------------------------------------------//
vector<bnsgaii::population> bnsgaii::fast_non_dominated(vector<population>& new_b_sol,vector<vector<double>> combine_sol){
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



// bnsgaii::solution bnsgaii::crowding_dis_assign(solution this_rank){
   
//     int solution_num=this_rank.size();
//     vector<double> init_dis(solution_num,0.0);
//     vector<vector<double>> this_obj_val(solution_num);
//     //population each_obj_val=fitness(F);
//     for(int i=0;i<this_rank.size();i++)
//         this_obj_val[i]=obj_val[this_rank[i]];

//     vector<int>tmp_id(solution_num);
//     vector<vector<int>> sorting_id(obj_num);
//     for(int i=0;i<solution_num;i++)
//         tmp_id[i]=i;
//     for(int i=0;i<obj_num;i++)
//         sorting_id[i]=tmp_id;
//     //列行互換
//     vector<vector<double>> row_obj_val(obj_num,vector<double>(solution_num));
//     for(int i=0;i<obj_num;i++)
//         for(int j=0;j<solution_num;j++)
//             row_obj_val[i][j]=this_obj_val[j][i];
//     //對各目標進行排序
//     //for(int i=0;i<obj_num;i++)
//     obj_sorting(row_obj_val[0],sorting_id[0]);
//     sorting_id[1]=sorting_id[0];
//     reverse(sorting_id[1].begin(),sorting_id[1].end());
//     vector<double> origin_obj1;
//     origin_obj1=row_obj_val[1];
//     for(int i=0;i<solution_num;i++)
//         row_obj_val[1][i]=origin_obj1[sorting_id[1][i]];
    
//     //記得要刪除-------------------------------------//
//     /*
//     cout<<"Have been sorted:"<<endl;
//     for(int i=0;i<obj_num;i++){
//         for(int j=0;j<solution_num;j++)
//             cout<<sorting_id[i][j]<<" ";
//         cout<<endl;
//     }*/
//     //---------------------------------------------------//
//     //記得要刪除-------------------------------------//
//     /*
//     for(int i=0;i<obj_num;i++){
//         for(int j=0;j<solution_num;j++)
//             cout<<row_obj_val[i][j]<<" ";
//         cout<<endl;
//     }*/
//     //---------------------------------------------------//
//     //記得要刪除-------------------------------------//
//     /*
//     cout<<"obj1:";
//     for(int i=0;i<solution_num;i++)
//         cout<<row_obj_val[0][i]<<" ";
//     cout<<endl;
//      cout<<"obj2:";
//     for(int i=0;i<solution_num;i++)
//         cout<<row_obj_val[1][i]<<" ";
//     cout<<endl<<endl;*/
//     //---------------------------------------------------//

//     //計算擁擠程度
//     init_dis[sorting_id[0][0]]=numeric_limits<double>::max();
//     init_dis[sorting_id[0][solution_num-1]]=numeric_limits<double>::max();
//     for(int i=0;i<obj_num;i++)
//         for(int j=1;j<solution_num-1;j++)
//             init_dis[sorting_id[0][j]]+=fabs(row_obj_val[i][j+1]-row_obj_val[i][j-1])/fabs(row_obj_val[i][solution_num-1]-row_obj_val[i][0]);
    
//     obj_sorting(init_dis,this_rank);
//     return this_rank;
// }
//-------------------------------------------Insertion sort--------------------------------------------------//
// void bnsgaii::obj_sorting(vector<double> &sort_obj,solution& sorting_id){
//     double max,tmp;
//     int max_pos;
//     for(int i=0;i<sort_obj.size();i++){
//         max=sort_obj[i];
//         max_pos=i;
//         for(int j=i+1;j<sort_obj.size();j++){
//             if(max<sort_obj[j]){
//                 max=sort_obj[j];
//                 max_pos=j;
//             }
//         }
//         //目標值排序
//         tmp=sort_obj[i];
//         sort_obj[i]=sort_obj[max_pos];
//         sort_obj[max_pos]=tmp;
//         //id排序
//         tmp=sorting_id[i];
//         sorting_id[i]=sorting_id[max_pos];
//         sorting_id[max_pos]=tmp;
//     }

// }
//------------------------------------------------------------------------------------------------------------//
vector<vector<double>> bnsgaii::fitness(vector<vector<double>> & all_sol) {
    if(problem_func=="SCH")
        return SCH(all_sol);
    else if(problem_func=="FON")
        return FON(all_sol);
    else if(problem_func=="POL")
        return POL(all_sol);
    else if(problem_func=="KUR")
        return KUR(all_sol);
    else if(problem_func=="ZDT1")
        return ZDT1(all_sol);
    else if(problem_func=="ZDT2")
        return ZDT2(all_sol);
    else if(problem_func=="ZDT3")
        return ZDT3(all_sol);
    else if(problem_func=="ZDT4")
        return ZDT4(all_sol);
    else if(problem_func=="ZDT6")
        return ZDT6(all_sol);
    else if(problem_func=="UF1")
        return UF1(all_sol);

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

bnsgaii::solution bnsgaii::crowding_dis_assign(solution this_rank){
   
    int front_size=this_rank.size();
    vector<double> init_dis(front_size,0.0);

    population sorting_id(obj_num,solution(front_size));
    population obj_id(obj_num);
    for(int i=0;i<obj_num;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort(i,obj_id[i],sorting_id[i],0,front_size-1);
    }
    
    //計算擁擠程度
    for (int i=0; i<obj_num; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
    
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<front_size-1;j++){
            if(init_dis[j]!=numeric_limits<double>::infinity()){
                if(obj_val[this_rank[sorting_id[i][front_size-1]]][i]== obj_val[this_rank[sorting_id[i][0]]][i])
                    init_dis[sorting_id[i][j]]+=0.0;
                else
                    init_dis[sorting_id[i][j]]+=(obj_val[this_rank[sorting_id[i][j+1]]][i]-obj_val[this_rank[sorting_id[i][j-1]]][i])/(obj_val[this_rank[sorting_id[i][front_size-1]]][i]-obj_val[this_rank[sorting_id[i][0]]][i]);
            }
        }
    //距離密度值進行正規化
    for (int j=0; j<front_size; j++)
        if (init_dis[j]!= numeric_limits<double>::infinity())
           init_dis[j] = init_dis[j]/obj_num;
        
    
    quick_sort_dis(init_dis,this_rank,0,front_size-1);
    reverse(this_rank.begin(),this_rank.end());
    return this_rank;
}
void bnsgaii::quick_sort_dis(vector<double>& init_dis,solution& this_rank,int left,int right){
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        temp = init_dis[right];
        init_dis[right] = init_dis[index];
        init_dis[index] = temp;

        temp = this_rank[right];
        this_rank[right] = this_rank[index];
        this_rank[index] = temp;

        pivot = init_dis[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (init_dis[j] <= pivot)
            {
                i+=1;
                temp = init_dis[j];
                init_dis[j] = init_dis[i];
                init_dis[i] = temp;
                
                temp = this_rank[j];
                this_rank[j] = this_rank[i];
                this_rank[i] = temp;
            }
        }
        index=i+1;
        temp = init_dis[index];
        init_dis[index] = init_dis[right];
        init_dis[right] = temp;

        temp = this_rank[index];
        this_rank[index] = this_rank[right];
        this_rank[right] = temp;

        quick_sort_dis(init_dis,this_rank, left, index-1);
        quick_sort_dis(init_dis,this_rank, index+1, right);
    }
}

#endif
       