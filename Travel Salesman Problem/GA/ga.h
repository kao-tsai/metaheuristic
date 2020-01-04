#ifndef __GA_H_INCLUDED__
#define __GA_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include<iomanip>
using namespace std;

class ga{
	public:
		typedef vector<int> solution;
        typedef vector<solution> population;
		ga(int,int,string,int,double,double,string,string);
		solution run();
	private:
		void init();
        double dis_cal(double,double,double,double);
        void generate_new_solutions();
		double evaluation(solution&);

        population PMX(population);
        population CX(population);
        population OX(population);

        void mutation(population&);
        vector<double> fitness(population);
		population RWS(population);
        population TS(population);

private:
	int numRuns;
	int numIter;
    int city_num;
	string filename;
    string selection;
    string cross_method;
    population node;
    vector<vector<double>> node_dis;
	population currentSol;
	vector<double> all_obj_value;
	int best_obj_value;
	solution best_sol;
    int population_num;
	double mutation_pro;
    double crossover_pro;
};

ga::ga(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
double  xCrossoverPro,
double  xMutationPro,
string xSelection,
string xCrossMethod
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	filename=xfilename;
    population_num=xPopulationNum;
    crossover_pro=xCrossoverPro;
    mutation_pro=xMutationPro;
    selection=xSelection;
    cross_method=xCrossMethod;
}
double ga::dis_cal(double x1,double y1,double x2,double y2){
	
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));

}

vector<double> ga::fitness(population all_sol) {
    vector<double> tmp(population_num);
	for (int i = 0; i < population_num; i++){
            tmp[i]=evaluation(all_sol[i]);
            if(best_obj_value>tmp[i]){
                    best_obj_value=tmp[i];
                    best_sol=all_sol[i];
            }
	}
    return tmp;
}
void ga::init(){
    string line,con;
    solution tmp;
    if (!filename.empty()){
        ifstream file(filename);
		node.clear();
		while(file>>line){
            file>>line;
            tmp.push_back(atoi(line.c_str()));
            file>>line;
            tmp.push_back(atoi(line.c_str()));
            node.push_back(tmp);
            tmp.clear();
        }
	file.close();
    city_num=node.size();
    node_dis=vector<vector<double>>(city_num);
    for(int i=0;i<city_num;i++)
			node_dis[i].resize(city_num);

	double result;
    for(int i=0;i<city_num;i++)
        for(int j=i;j<city_num;j++)
        {		
                if(i!=j){
                    result=dis_cal(node[i][0],node[i][1],node[j][0],node[j][1]);
                    node_dis[i][j]=result;
                    node_dis[j][i]=result;
                }
                else{
                    node_dis[i][j]=0;		
                }            
        }
    }
    else{
		cout<<"Need Input!"<<endl;
    }

}

void ga::generate_new_solutions(){
    best_obj_value=10000;
    currentSol.clear();
    currentSol.resize(population_num);
    solution city_set;
    solution tmp;
    int random_city;
    for(int i=0;i<city_num;i++)
        city_set.push_back(i);

    for (int i = 0; i < population_num; i++){
        tmp=city_set;
        for(int j=0;j<city_num;j++)
        {
            random_city=rand()%tmp.size();
            currentSol[i].push_back(tmp[random_city]);
            tmp.erase(tmp.begin()+random_city);
        }
    }
}

//----------------------------------------------------------------------Selection-------------------------------------------------------------------------------//
ga::population ga::RWS(population all_sol){
    double  do_selection,min;
    int min_position,prevision=100,sum=0;
    vector<double> wheel(population_num);
    vector<double> restore(population_num);
    population new_pop;
    restore=all_obj_value;
    //設置Rank
    for(int i=0;i<population_num;i++)
    {
        
        min=all_obj_value[0];
        min_position=0;
        for(int j=1;j<population_num;j++)
        {
            if(min>all_obj_value[j])
            {
                min=all_obj_value[j];
                min_position=j;
            }
            //sum+=all_obj_value[i];
        }
        wheel[min_position]=i+1;
        sum+=i+1;
        all_obj_value[min_position]=10000;
    }
    
    all_obj_value=restore;
    cout<<"這是sum:"<<sum<<endl;
    //計算輪盤機率
     for(int i=0;i<population_num;i++){
        if(i!=0)
            wheel[i]=wheel[i]+wheel[i-1];
        //cout<<wheel[i]<<endl;
     }
     /*
    for(int i=0;i<population_num;i++){
        wheel[i]=wheel[i]*prevision;
        cout<<wheel[i]<<endl;
     }
    */
    for(int i=0;i<population_num;i++)
    {
            do_selection=rand()%sum;
            //do_selection=rand()%prevision;
            for(int j=0;j<population_num;j++)
                if(do_selection<wheel[j]){
                    new_pop.push_back(all_sol[j]);   
                    break;
                }
    }
    
    return new_pop;
}

ga::population ga::TS(population all_sol){
    population new_pop;
    bool flag;
    int tournament_size=3,tmp,candidate_best;
    solution candidate;
    double min;
    for(int i=0;i<population_num;i++){
        candidate.clear();
        do{
        flag=true;
        tmp=rand()%population_num;

        for(int j=0;j<candidate.size();j++)
            if(tmp==candidate[j]){
                flag=~flag;
                break;
            }
        if(flag)
            candidate.push_back(tmp);  
        }while(candidate.size()<tournament_size);
        
        min=evaluation(all_sol[candidate[0]]);
        
        candidate_best=candidate[0];
        for(int j=1;j<tournament_size;j++){
            tmp=all_obj_value[candidate[j]];
            if(tmp<min)
            {
                min=tmp;
                candidate_best=candidate[j];
            }

        }
        new_pop.push_back(all_sol[candidate_best]);
    }
    return new_pop;
}

//----------------------------------------------------------------------Crossover-------------------------------------------------------------------------------//
ga::population ga::PMX(population all_sol){
    int cut1,cut2,parent1,parent2,tmp,change1,change2;
    solution temp_sol(city_num);
    solution temp_list;
    population list;
    population new_pop;

    while(all_sol.size()>0){
            
        if(all_sol.size()==1){
            new_pop.push_back(all_sol[0]);
            break;
        }

        parent1=rand()%all_sol.size();
        do{
            parent2=rand()%all_sol.size();
        }while(parent2==parent1);

        temp_sol=all_sol[parent1];


        //crossover probability------------------------------------------------------------------------------------------------------------------------------------------------
        double r=(double)rand()/RAND_MAX;
        
        if((r<crossover_pro) && (all_sol[parent1]!=all_sol[parent2])){
           
            //選擇切點------------------------------------------------------------------------------------------------
            cut1=rand()%(city_num);
            do{
                cut2=rand()%(city_num);
            }while(cut1==cut2);

            if(cut1>cut2){
                tmp=cut1;
                cut1=cut2;
                cut2=tmp;
            }
        
            //建立list--------------------------------------------------------------------------------------------------
            list.clear();
            for (int x = cut1; x <= cut2; x++){
                change1=-1;
                change2=-1;
                temp_list.clear();
                temp_list.push_back(all_sol[parent1][x]);
                temp_list.push_back(all_sol[parent2][x]);
                if(temp_list[0]==temp_list[1])
                {
                    temp_list.clear();
                    continue;
                }
                //整理List
                for(int k=0;k<list.size();k++){
                    //檢查頭
                    if(temp_list[0]==list[k][1])
                    {
                        list[k][1]=temp_list[1];
                        if(list[k][0]==list[k][1])
                        {
                            list.erase(list.begin()+k);
                            change1=-2;
                            break;
                        }
                        change1=k;
                    }
                    //檢查尾
                    if(temp_list[1]==list[k][0])
                    {
                        list[k][0]=temp_list[0];
                        if(list[k][0]==list[k][1])
                        {
                            list.erase(list.begin()+k);
                            change2=-2;
                            break;
                        }
                        change2=k;
                    }
                }
                if(change1!=-1 && change2!=-1) {
                        list[change1][1]=list[change2][1];
                        list.erase(list.begin()+change2);
                }
                else if(change1==-1 && change2==-1)
                    list.push_back(temp_list);      
                temp_list.clear();
            }
            //印出list
            // cout<<list.size()<<endl;
            // for(int m=0;m<list.size();m++){
            //     for(int n=0;n<list[m].size();n++)
            //         cout<<list[m][n]<<" ";
            //     cout<<endl;
            // }
            // cout<<"this"<<endl;
            //做crossover---------------------------------------------------------------------------------------------------
            for (int j = cut1; j <= cut2; j++){
                    all_sol[parent1][j] = all_sol[parent2][j];
                    all_sol[parent2][j] = temp_sol[j];
            }
            //染色體修正cut1前-------------------------------------------------------------------------------------------
            
            for(int i=0;i<cut1;i++)
            {
                //尋找修正基因
                for(int j=0;j<list.size();j++)
                {
                    if(list[j][0]==all_sol[parent1][i])
                        all_sol[parent1][i]=list[j][1];
                       
                    else if(list[j][1]==all_sol[parent1][i])
                        all_sol[parent1][i]=list[j][0];
                       
                    if(list[j][0]==all_sol[parent2][i])
                       all_sol[parent2][i]=list[j][1];
                       
                    
                    else if(list[j][1]==all_sol[parent2][i])
                        all_sol[parent2][i]=list[j][0];
                  
                }
            }        
            //染色體修正cut2後------------------------------------------------------------------------------------------
            for(int i=cut2+1;i<city_num;i++)
            {
                //尋找修正基因
                for(int j=0;j<list.size();j++)
                {
                    if(list[j][0]==all_sol[parent1][i])
                        all_sol[parent1][i]=list[j][1];
                        
                    
                    else if(list[j][1]==all_sol[parent1][i])
                    
                        all_sol[parent1][i]=list[j][0];
                    
                    

                    if(list[j][0]==all_sol[parent2][i])
                    
                        all_sol[parent2][i]=list[j][1];
                     
                    
                    else if(list[j][1]==all_sol[parent2][i])
                    
                        all_sol[parent2][i]=list[j][0];
                           
                }
                
            }
            
        }
        
        new_pop.push_back(all_sol[parent1]);
        new_pop.push_back(all_sol[parent2]);
        
        all_sol.erase(all_sol.begin()+parent1);
        all_sol.erase(all_sol.begin()+parent2);
    }
    return new_pop;
}

ga::population ga::CX(population all_sol){
    int parent1,parent2,find_next,start_position;
    solution need_crossover;
    solution temp_sol;
    population new_pop;

    while(all_sol.size()>0){
        need_crossover=solution(city_num,1);
        if(all_sol.size()==1){
            new_pop.push_back(all_sol[0]);
            break;
        }
        parent1=rand()%all_sol.size();
        do{
            parent2=rand()%all_sol.size();
        }while(parent2==parent1);
        temp_sol=all_sol[parent1];
        double r=(double)rand()/RAND_MAX;
        if(r<crossover_pro){
            
            start_position=rand()%city_num;
            find_next=start_position;
            need_crossover[start_position]=0;
            //排除不需要crossover的基因（建立cycle）
            while(all_sol[parent1][start_position]!=all_sol[parent2][find_next])
                for(int i=0;i<city_num;i++)
                    if(all_sol[parent1][i]==all_sol[parent2][find_next])
                    {
                        need_crossover[i]=0;
                        find_next=i;
                        break;
                    }

            
            //做crossover（選0的部份）
            for(int i=0;i<city_num;i++)
                if(need_crossover[i])
                {
                    all_sol[parent1][i]=all_sol[parent2][i];
                    all_sol[parent2][i]=temp_sol[i];
                }
                
        }
        new_pop.push_back(all_sol[parent1]);
        new_pop.push_back(all_sol[parent2]);
        all_sol.erase(all_sol.begin()+parent1);
        all_sol.erase(all_sol.begin()+parent2);
    }
    return new_pop;
}
ga::population ga::OX(population all_sol){
    int cut1,cut2,parent1,parent2,tmp;
    solution p1_cross_gene;
    solution p2_cross_gene;
    population new_pop;
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    while(all_sol.size()>0){
            
        if(all_sol.size()==1){
            new_pop.push_back(all_sol[0]);
            break;
        }

        parent1=rand()%all_sol.size();
        do{
            parent2=rand()%all_sol.size();
        }while(parent2==parent1);

        p1_cross_gene=all_sol[parent2];
        p2_cross_gene=all_sol[parent1];

        //crossover probability------------------------------------------------------------------------------------------------------------------------------------------------

        double r=(double)rand()/RAND_MAX;
        if(r<crossover_pro){

            //選擇切點------------------------------------------------------------------------------------------------
            cut1=rand()%(city_num);
            do{
                cut2=rand()%(city_num);
            }while(cut1==cut2);

            if(cut1>cut2){
                tmp=cut1;
                cut1=cut2;
                cut2=tmp;
            }
            //清除重複的基因--------------------------------------------------------------------------
            for(int i=cut1;i<=cut2;i++)
            {
                for(int j=0;j<p1_cross_gene.size();j++)
                    if(p1_cross_gene[j]==all_sol[parent1][i]){
                        p1_cross_gene.erase(p1_cross_gene.begin()+j);  
                        break;
                    }
                for(int j=0;j<p2_cross_gene.size();j++)
                    if(p2_cross_gene[j]==all_sol[parent2][i]) {
                        p2_cross_gene.erase(p2_cross_gene.begin()+j);
                        break;
                    }
            }

            //填入基因-----------------------------------------------------------------------------------
           for(int i=cut2;i>=cut1;i--){
                p1_cross_gene.insert(p1_cross_gene.begin()+cut1,all_sol[parent1][i]);
                p2_cross_gene.insert(p2_cross_gene.begin()+cut1,all_sol[parent2][i]);
           }
        }
        new_pop.push_back(p1_cross_gene);
        new_pop.push_back(p2_cross_gene);
        all_sol.erase(all_sol.begin()+parent1);
        all_sol.erase(all_sol.begin()+parent2);
    }
    return new_pop;
}
//----------------------------------------------------------------------Mutation--------------------------------------------------------------------------------//
void ga::mutation(population &all_sol){
    int mutation_city1,mutation_city2,tmp;
    
    for(int i=0;i<population_num;i++){
            double r=(double)rand()/RAND_MAX;
            if(r<mutation_pro){
                mutation_city1=rand()%city_num;
                do{
                    mutation_city2=rand()%city_num;
                }while(mutation_city1==mutation_city2);

                tmp=all_sol[i][mutation_city1];
                all_sol[i][mutation_city1]=all_sol[i][mutation_city2];
                all_sol[i][mutation_city2]=tmp;
                
            }
    }
}

ga::solution ga::run(){

    double avg_obj_value;
    vector<double>iter_obj_avg(numIter,0.0);
    init();
    
    for(int i=0;i<numRuns;i++){

        generate_new_solutions();
        
        avg_obj_value=0.0;

        for(int j=0;j<numIter;j++){

            all_obj_value=fitness(currentSol);
            
            if(selection=="ts")
                currentSol=TS(currentSol);
            else if(selection=="rws")
                currentSol=RWS(currentSol);
            
            if(cross_method=="pmx")
                currentSol=PMX(currentSol);
            else if(cross_method=="cx")
                currentSol=CX(currentSol);
            else if(cross_method=="ox")
                currentSol=OX(currentSol);

            
            mutation(currentSol);
          
            
            iter_obj_avg[j]+=best_obj_value;
            
        }
        
    }
    
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;
    
    return best_sol;

}

double ga::evaluation(solution& sol){

    double  distance=0;
    for(int i=0;i<city_num-1;i++)
        distance+=node_dis[sol[i]][sol[i+1]];
    distance+=node_dis[sol[city_num-1]][sol[0]];
    
    return distance;
}
#endif
       