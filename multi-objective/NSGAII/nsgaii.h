#ifndef __NSGAII_H_INCLUDED__
#define __NSGAII_H_INCLUDED__
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

class nsgaii{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		nsgaii(int,int,int,double,double,string);
		population run();
	private:
		void init(population&);
		void crossover(population&);
        void mutation(population&);
        population fitness(population&);
		population RWS(population);
        void TS(population&);
        population fast_non_dominated(population);
        vector<int> crowding_dis_assign(vector<int>);
        void obj_sorting(solution&,vector<int>&);

       

private:
	int numRuns;
	int numIter;
	population obj_val;
	solution best_obj_value;
    int population_num;
	double mutation_pro;
    double crossover_pro;
    string problem_func;
    int dimension; //每一個解的維度
    int obj_num; //目標數量
    solution upperbound;
    solution lowerbound;
    double dimm;
    double dimc;
    vector<int> nbits;
    
};
nsgaii::nsgaii(int xNumRuns,
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
    population_num=xPopulationNum;
    problem_func=xproblem_func;
    func_init(problem_func);
    dimension=idimension;
    obj_num=iobj_num;
    upperbound=ub;
    lowerbound=lb;
    dimm=20;
    dimc=20;
}


void nsgaii::init(population& currentSol){
    //best_obj_value=0;
	for (int i = 0; i < population_num; i++)
        for(int j=0;j<dimension;j++)
            currentSol[i][j]=(upperbound[j]-lowerbound[j])*((double)rand()/RAND_MAX)+lowerbound[j];
}
/*
nsgaii::population nsgaii::RWS(population all_sol){
    int sum=0;
    vector<double> wheel(population_num);
    population new_pop;
    int selection;
    for(int i=0;i<population_num;i++)
    {
        sum+=obj_value[i];
        wheel[i]=sum;
    }

    for(int i=0;i<population_num;i++)
    {
            selection=rand()%sum;
            for(int j=0;j<population_num;j++)
                if(selection<wheel[i])
                    new_pop.push_back(all_sol[j]);
    }
    return new_pop;
}
*/

void nsgaii::TS(population &all_sol){
    population new_pop(population_num);

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

// void nsgaii::crossover(population& all_sol){
//     int y;
//     double rnd,par1,par2,chld1,chld2,betaq,beta,alpha;
//     double y1,y2,yu,yl;
//     y=0;
    
//     for(int i = 0; i < population_num/2; i++){
//         rnd =(double)rand()/RAND_MAX;
      
//         /*Check Whether the cross-over to be performed*/
//         if(rnd <= crossover_pro){
           
// 	        /*Loop over no of variables*/
// 	        for(int j = 0;j < dimension;j++){
//             /*Selected Two Parents*/  
// 	            rnd = (double)rand()/RAND_MAX;             
//                 /* Check whether variable is selected or not*/
//                 if(rnd <= 0.5){
//                     /*Variable selected*/
//                     if(fabs(par1 - par2) > EPS){
//                         if(par2 > par1){
//                             y2 = par2;
//                             y1 = par1;
//                         }
//                         else{
//                             y2 = par1;
//                             y1 = par2;
//                         }
//                         yl = lowerbound[j];
//                         yu = upperbound[j];
                       
//                         rnd=(double)rand()/RAND_MAX;
//                         beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
//                         alpha = 2.0 - pow(beta,-(dimc+1.0));
                        
//                         if (rnd <= 1.0/alpha)
//                         {
//                            betaq = pow ((rand*alpha),(1.0/(dimc+1.0)));
//                         }
//                         else
//                         {
//                             betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
//                         }
//                         c1 = 0.5*((y1+y2)-betaq*(y2-y1));
//                         beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
//                         alpha = 2.0 - pow(beta,-(eta_c+1.0));
//                          if (rand <= (1.0/alpha))
//                         {
//                             betaq = pow ((rand*alpha),(1.0/(dimc+1.0)));
//                         }
//                         else
//                         {
//                             betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(dimc+1.0)));
//                         }
//                          c2 = 0.5*((y1+y2)+betaq*(y2-y1));
//                         if (c1<yl)
//                             c1=yl;
//                         if (c2<yl)
//                             c2=yl;
//                         if (c1>yu)
//                             c1=yu;
//                         if (c2>yu)
//                             c2=yu;
//                         rnd=(double)rand()/RAND_MAX;
//                         if (rnd<=0.5)
//                         {
//                             child1->xreal[i] = c2;
//                             child2->xreal[i] = c1;
//                         }
//                         else
//                         {
//                             child1->xreal[i] = c1;
//                             child2->xreal[i] = c2;
//                         }
// 		        }
//                 else{
//                     /*Copying the children to parents*/
//                     chld1 = par1;
//                     chld2 = par2;
//                 }
//                 all_sol[y][j]= chld1;
//                 all_sol[y+1][j] = chld2;
// 	        }
// 	    }
//         y=y+2;
//     }
  
// }

void nsgaii::crossover(population& all_sol){

}
/*
void nsgaii::mutation(population& all_sol){
    double rnd,delta,indi,deltaq;
    double y,yl,yu,val,xy;
    for(int i=0;i<population_num;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            rnd=(double)rand()/RAND_MAX;
            if(rnd<=mutation_pro)
            {
                y=all_sol[i][j];
                yl=lowerbound[j];
                yu=upperbound[j];
                //計算delta
                if(y>yl)
                {
                    if((y-yl) < (yu-y))
                        delta = (y - yl)/(yu - yl);
                    else
                        delta = (yu - y)/(yu-yl);

                    rnd=(double)rand()/RAND_MAX;
                    indi=1.0/(dimm+1.0);
                    if(rnd<=0.5)
                    {
                        xy=1.0-delta;
                        val=2*rnd+(1-2*rnd)*(pow(xy,(dimm+1)));
                        deltaq=pow(val,indi)-1.0;
                    }
                    else
                    {
                        xy=1.0-delta;
                        val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(dimm+1)));
		                deltaq = 1.0 - (pow(val,indi));
                    }
                    y = y + deltaq * (yu-yl);
                    if (y < yl) y=yl; 
                    if (y > yu) y=yu;
                    all_sol[i][j]=y;
                }
                else
                {
                    xy=(double)rand()/RAND_MAX;
                    all_sol[i][j]=xy*(yu-yl)+yl;
                }             
            }
        }
    }
}*/

//-----------------------------------------------------run nsga2 algorithm----------------------------------------------------------//
nsgaii::population nsgaii::run(){
    
    vector<double>iter_obj_avg(numIter,0.0);
    population currentSol(population_num,solution(dimension));
    population newSol;
   
    for(int i=0;i<numRuns;i++){
        
        init(currentSol);
        
        //印出所有人口測試
        /*
        for(int a=0;a<population_num;a++){
            for(int b=0;b<dimension;b++)
                cout<<currentSol[a][b]<<" ";
            cout<<endl;
        }*/
        
       cout<<"this is currentSol size:"<<currentSol[0].size()<<endl;
        newSol=fast_non_dominated(currentSol);
        cout<<"this is newSol size:"<<newSol[0].size()<<endl;
        for(int j=0;j<numIter;j++){
            
            TS(newSol);
            crossover(newSol);
            
            mutation(newSol);
            
            //IGD
            newSol.insert(newSol.end(),currentSol.begin(),currentSol.end());//combination
            
            currentSol=fast_non_dominated(newSol);
            
            newSol=currentSol;
            //iter_obj_avg[j]+=best_obj_value;
        }
        obj_val=fitness(newSol);
    }
    /*
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
    */
    return obj_val;

}
//------------------------------------------------------------------------------------------------------------------------------------------//
nsgaii::population nsgaii::fast_non_dominated(population combine_sol){
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
        
            if(obj_val[i][0]<=obj_val[j][0]&&obj_val[i][1]<=obj_val[j][1])
                if(obj_val[i][0]==obj_val[j][0]&&obj_val[i][1]==obj_val[j][1])
                    continue;
                else
                {
                    set[i].push_back(j);
                    set[j][0]++;
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
    vector<vector<int>> rank;
    vector<int> temp_rank;
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
    }
    */
    //-----------------------------------//

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
        rank[rank.size()-1].pop_back();
        new_population--;
    }
     
    //最終排序人口回傳
    population final_population(population_num);
    int x=0;
    for(int i=0;i<rank.size();i++){
        for(int j=0;j<rank[i].size();j++){
            final_population[x]=combine_sol[rank[i][j]];
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



vector<int> nsgaii::crowding_dis_assign(vector<int> this_rank){
    int solution_num=this_rank.size();
    solution init_dis(solution_num,0.0);
    population this_obj_val(solution_num);
    //population each_obj_val=fitness(F);
    for(int i=0;i<this_rank.size();i++)
        this_obj_val[i]=obj_val[this_rank[i]];

    vector<int>tmp_id(solution_num);
    vector<vector<int>> sorting_id(obj_num);
    for(int i=0;i<solution_num;i++)
        tmp_id[i]=i;
    for(int i=0;i<obj_num;i++)
        sorting_id[i]=tmp_id;
    //列行互換
    population row_obj_val(obj_num,solution(solution_num));
    for(int i=0;i<obj_num;i++)
        for(int j=0;j<solution_num;j++)
            row_obj_val[i][j]=this_obj_val[j][i];
    //對各目標進行排序
    for(int i=0;i<obj_num;i++)
        obj_sorting(row_obj_val[i],sorting_id[i]);
    
    //計算擁擠程度
    init_dis[sorting_id[0][0]]=numeric_limits<double>::max();
    init_dis[sorting_id[0][solution_num-1]]=numeric_limits<double>::max();
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<solution_num-1;j++)
            init_dis[sorting_id[i][j]]+=fabs(row_obj_val[i][j+1]-row_obj_val[i][j-1])/(row_obj_val[i][solution_num-1]-row_obj_val[i][0]);
    obj_sorting(init_dis,this_rank);
    return this_rank;
}
void nsgaii::obj_sorting(solution &sort_obj,vector<int>& sorting_id){
    double max,tmp;
    int max_pos;
    for(int i=0;i<sort_obj.size();i++){
        max=sort_obj[i];
        max_pos=i;
        for(int j=i+1;j<sort_obj.size();j++){
            if(max<sort_obj[j]){
                max=sort_obj[j];
                max_pos=j;
            }
        }
        //目標值排序
        tmp=sort_obj[i];
        sort_obj[i]=sort_obj[max_pos];
        sort_obj[max_pos]=tmp;
        //id排序
        tmp=sorting_id[i];
        sorting_id[i]=sorting_id[max_pos];
        sorting_id[max_pos]=tmp;
    }

}
nsgaii::population nsgaii::fitness(population& all_sol) {
    if(problem_func=="SCH")
        return SCH(all_sol);
    else if(problem_func=="FON")
        return FON(all_sol);
    else if(problem_func=="POL")
        return POL(all_sol);
}




#endif
       