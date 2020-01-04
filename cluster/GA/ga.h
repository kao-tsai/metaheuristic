#ifndef __GA_H_INCLUDED__
#define __GA_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include<iomanip>
using namespace std;

class ga{
	public:
		typedef vector<double> solution;
        typedef vector<solution> centroid;
        typedef vector<centroid> population;
		ga(int,int,string,int,int);
		centroid run();
	private:
		void init();
        double dis_cal(solution,solution);
        void generate_new_solutions();
		double evaluation(centroid&);
        double frand(double,double);
        void mutation(population&);
        void fitness(population&);
		population RWS(population);
        population TS(population);
        population crossover(population);
private:
	int numRuns;
	int numIter;
    int dim;
	string filename;
    centroid  node;
	population currentSol;
	solution all_obj_value;
	int best_obj_value;
	centroid best_sol;
    int population_num;
    int cluster_num;
    double dmax;
    double dmin;
    double crossover_pro;
    double mutation_pro;
};

ga::ga(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
int x_cluster_num
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	filename=xfilename;
    population_num=xPopulationNum;
    cluster_num=x_cluster_num;
}

void ga::fitness(population& sol) {
    double tmp_best;
	for (int i = 0; i < population_num; i++){
            tmp_best=evaluation(sol[i]);
            if(tmp_best<best_obj_value){
                    best_obj_value=tmp_best;
                    best_sol=sol[i];
            }
	}
}
/*
void ga::init(){
    string line,con;
    solution tmp;
    int count_dim=0;
    if (!filename.empty()){
        ifstream file("GA/"+filename);
		node.clear();
        stringstream ss;
        //決定問題的維度
        getline(file,line);
        
        ss<<line;
        while(getline(ss,con,',')){    
            count_dim++;
        }
        dim=count_dim-1;
        //紀錄各維度最大最小值
        ss.clear();
        ss=stringstream(line);
        dim_minmax.resize(dim);
        //abalone
        // getline(ss,con,',');
        //abalone
        for(int i=0;i<dim;i++){
            getline(ss,con,',');
            tmp.push_back((double)atof(con.c_str()));
            //存最小值
            dim_minmax[i].push_back(tmp[i]);
            //存最大值
            dim_minmax[i].push_back(tmp[i]);
        }
        node.push_back(tmp);
        
        //建立初始節點
        
		while(getline(file,line)){
            ss.clear();
            ss=stringstream(line);
            //abalone
            //getline(ss,con,',');
            //abalone
            for(int i=0;i<dim;i++){   
                getline(ss,con,',');
                tmp[i]=(double)atof(con.c_str());
            }
            node.push_back(tmp);
            //開始尋找各維度最大最小值----------
            for(int i=0;i<dim;i++){
                //存最小值
                if(tmp[i]<dim_minmax[i][0])
                    dim_minmax[i][0]=tmp[i];
                //存最大值
                else if(tmp[i]>dim_minmax[i][1])
                    dim_minmax[i][1]=tmp[i];
            }
        }
	file.close();
  
    //currentSol計得初始化;
    currentSol.resize(population_num);
    for(int i=0;i<population_num;i++){
        currentSol[i].resize(cluster_num);
        for(int j=0;j<cluster_num;j++)
            currentSol[i][j].resize(dim);
    }
    
    }
    else{
		cout<<"Need Input!"<<endl;
    }
    
}
*/
void ga::init(){
    string line,con;
    solution tmp;
    int count_dim=0;
    if (!filename.empty()){
        ifstream file("GA/"+filename);
		node.clear();
        stringstream ss;
        //決定問題的維度
        getline(file,line);
        
        ss<<line;
        while(getline(ss,con,','))
            count_dim++;
        dim=count_dim-1;
        //紀錄各維度最大最小值
        ss.clear();
        ss=stringstream(line);
     
        for(int i=0;i<dim;i++){
            getline(ss,con,',');
            tmp.push_back((double)atof(con.c_str()));
            if(i==0){
                dmin=tmp[i];
                dmax=tmp[i];
            }
           else
                if(tmp[i]<dmin)
                    dmin=tmp[i];
                else if(tmp[i]>dmax)
                    dmax=tmp[i]; 
        }
        node.push_back(tmp);
        
        //建立初始節點
		while(getline(file,line)){
            ss.clear();
            ss=stringstream(line);
            for(int i=0;i<dim;i++){   
                getline(ss,con,',');
                tmp[i]=(double)atof(con.c_str());
            }
            node.push_back(tmp);
            //開始尋找各維度最大最小值----------
            for(int i=0;i<dim;i++){
                //存最小值
                if(tmp[i]<dmin)
                    dmin=tmp[i];
                //存最大值
                else if(tmp[i]>dmax)
                   dmax=tmp[i];
            }
        }
        file.close();
    
        //currentSol計得初始化;
        currentSol.resize(population_num);
        for(int i=0;i<population_num;i++){
            currentSol[i].resize(cluster_num);
            for(int j=0;j<dim;j++)
                currentSol[i][j].resize(dim);
        }

    }
    else{
		cout<<"Need Input!"<<endl;
    }
    
}

/*
void ga::generate_new_solutions(){
    
    for (int i = 0; i < population_num; i++){
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++){
                currentSol[i][j][k]=frand(dim_minmax[k][0],dim_minmax[k][1]);
            }
    }

    best_obj_value=evaluation(currentSol[0]);
    best_sol=currentSol[0];
}
*/

void ga::generate_new_solutions(){
    solution rand_centroid;
    solution tmp;
    for (int i = 0; i < population_num; i++){
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++)
                currentSol[i][j][k]=frand(dmin,dmax);
    }
    best_obj_value=evaluation(currentSol[0]);
    best_sol=currentSol[0];
}
double ga::evaluation(centroid &sol){
    double min;
    double sum=0.0;
    double tmp_dis;
    for(int i=0;i<node.size();i++){
        min=dis_cal(node[i],sol[0]);
        for(int j=1;j<cluster_num;j++){
            if(min>(tmp_dis=dis_cal(node[i],sol[j])))
                min=tmp_dis;
        }
        sum+=min;
    }
    return sum;
}
//計算兩點之間的距離
double ga::dis_cal(solution p1,solution p2){
    double sum=0.0;
    for(int i=0;i<dim;i++)
        sum+=pow(p1[i]-p2[i],2.0);
    
    return sum;
}

double ga::frand(double fmin,double fmax){
    double f = (double)rand() / RAND_MAX;
    return fmin + f * (fmax - fmin);
}


//----------------------------------------------------------------------Selection-------------------------------------------------------------------------------//
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
ga::population ga::crossover(population all_sol){
    int parent1,parent2,find_next,start_position;

    solution temp_sol;
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
//----------------------------------------------------------------------Mutation--------------------------------------------------------------------------------//


ga::centroid ga::run(){
    int count=0;
    double tmp_best_obj_value;
    double avg_obj_value;
    population tmp_sol;
    vector<double>iter_obj_avg(numIter,0.0);
    
    init();
   
    for(int i=0;i<numRuns;i++){
        count++;
        
        generate_new_solutions();


        avg_obj_value=0.0;

        for(int j=0;j<numIter;j++){


            currentSol=tmp_sol;
            
            fitness(currentSol);
            iter_obj_avg[j]+=best_obj_value;
            
        }
         
    }
    
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;
    
    return best_sol;

}


#endif
       