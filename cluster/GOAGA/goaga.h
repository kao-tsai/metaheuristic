#ifndef __GOAGA_H_INCLUDED__
#define __GOAGA_H_INCLUDED__
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

class goaga{
	public:
		typedef vector<double> solution;
        typedef vector<solution> centroid;
        typedef vector<centroid> population;
		goaga(int,int,string,int,int,double,double,double,double);
		centroid run();
	private:
		void init();
        double dis_cal(solution,solution);
        void generate_new_solutions();
		double evaluation(centroid);
        double frand();
        double deg_forces(double);
        centroid transition(population,int,double);
        double coef_c(int);

        void fitness(population&);
        population TS(population);
        population crossover(centroid);
        centroid vec_dim_transfer(population);
        population vec_dim_transfer(centroid);
        void mutation(population&);

private:
	int numRuns;
	int numIter;
    int dim;
	string filename;
    centroid  node;
	population currentSol;
	solution all_obj_value;
	double best_obj_value;
	centroid best_sol;
    int population_num;
    int cluster_num;

    double dmax;
    double dmin;
    double cmax;
    double cmin;
    double crossover_pro;
    double mutation_pro;
};

goaga::goaga(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
int x_cluster_num,
double xcrossover_pro,
double xmutation_pro,
double xcmax,
double xcmin
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	filename=xfilename;
    population_num=xPopulationNum;
    cluster_num=x_cluster_num;
    crossover_pro=xcrossover_pro;
    mutation_pro=xmutation_pro;
    cmax=xcmax;
    cmin=xcmin;
}

double goaga::deg_forces(double r){
    double f=0.5;
    double l=1.5;
    return f*exp((-r)/l)-exp(-r);
}

void goaga::fitness(population &sol) {
    double tmp_best;
	for (int i = 0; i < population_num; i++){
            tmp_best=evaluation(sol[i]);
            all_obj_value[i]=tmp_best;
            if(tmp_best<best_obj_value){
                    best_obj_value=tmp_best;
                    best_sol=sol[i];
            }
	}
}


void goaga::init(){
    string line,con;
    solution tmp;
    int count_dim=0;
    if (!filename.empty()){
        ifstream file("dataset/"+filename);
		node.clear();
        stringstream ss;
        //決定問題的維度
        getline(file,line);
        ss<<line;
        while(getline(ss,con,','))
            count_dim++;
        dim=count_dim;
        //紀錄各維度最大最小值
        ss.clear();
        ss=stringstream(line);
        //abalone
        //getline(ss,con,',');
        //abalone
        for(int i=0;i<dim;i++){
            getline(ss,con,',');
            tmp.push_back((double)stod(con));
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
            // abalone
            //getline(ss,con,',');
            // abalone
            for(int i=0;i<dim;i++){   
                
                getline(ss,con,',');
                tmp[i]=(double)stod(con);
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
        currentSol=population(population_num,centroid(cluster_num,solution(dim)));
        all_obj_value=solution(population_num);
    }
    else{
		cout<<"Need Input!"<<endl;
    }
  
}


void goaga::generate_new_solutions(){
    solution rand_centroid;
    solution tmp;
    for (int i = 0; i < population_num; i++)
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++)
                currentSol[i][j][k]=frand();
    
    best_obj_value=evaluation(currentSol[0]);
    best_sol=currentSol[0];
}
double goaga::evaluation(centroid sol){
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
double goaga::dis_cal(solution p1,solution p2){
    double sum=0.0;
    for(int i=0;i<dim;i++)
        sum+=pow(p1[i]-p2[i],2.0);
    return sum;
}

double goaga::frand(){
    double f = (double)rand()/ RAND_MAX;
    return dmin + f * (dmax - dmin);
}
double goaga::coef_c(int cur_iter){
    double a;
    a=cmax-((cmax-cmin)/(double)(numIter-1))*(double)cur_iter;
    return a;
}

//----------------------------------------------------------------------force-------------------------------------------------------------------------------------//


goaga::centroid goaga::transition(population sol,int sol_num,double c){

    double normalize_val;
    double tmp_dis;
    centroid tmp_sol=centroid(cluster_num,solution(dim));

    for(int i=0;i<population_num;i++){
        if(sol_num==i)
            continue;

        for(int j=0;j<cluster_num;j++){
            tmp_dis=sqrt(dis_cal(sol[sol_num][j],sol[i][j]));    
            // if(tmp_dis==0)
            //     tmp_dis=0.00001;
            for(int k=0;k<dim;k++){
                //作正規化
                normalize_val=((fabs(sol[i][j][k]-sol[sol_num][j][k])-dmin)/(dmax-dmin))*3+1;
                tmp_sol[j][k]+=c*((dmax-dmin)/2.0)*deg_forces(normalize_val)*((sol[i][j][k]-sol[sol_num][j][k])/tmp_dis);
            }    
        }
    }
    
    //此處best_sol有爭議不確定用單點或是向量
    for(int i=0;i<cluster_num;i++)
        for(int j=0;j<dim;j++){

            tmp_sol[i][j]=c*tmp_sol[i][j]+best_sol[i][j];
            
            bool flag=true;
            while(flag){
                flag=false;
                if(tmp_sol[i][j]<dmin)
                {
                    tmp_sol[i][j]=tmp_sol[i][j]+(dmax-dmin)*((double)rand()/RAND_MAX);
                    flag=true;
                }
                else if(tmp_sol[i][j]>dmax)
                {
                    tmp_sol[i][j]= tmp_sol[i][j]-(dmax-dmin)*((double)rand()/RAND_MAX);
                    flag=true;
                }
            } 
        }

   return tmp_sol;
}


//----------------------------------------------------------------------Selection-------------------------------------------------------------------------------//
goaga::population goaga::TS(population all_sol){
    population new_pop;
    bool flag;
    int tournament_size=3,tmp,candidate_best;
    vector<int> candidate;
    double tmp_min;
    for(int i=0;i<population_num;i++){
        candidate.clear();
        do{
        flag=true;
        tmp=rand()%population_num;

        for(int j=0;j<candidate.size();j++)
            if(tmp==candidate[j]){
                flag=false;
                break;
            }
        if(flag)
            candidate.push_back(tmp);  
        }while(candidate.size()<tournament_size);
        
        tmp_min=all_obj_value[candidate[0]];
        candidate_best=candidate[0];

        for(int j=1;j<tournament_size;j++){
            
            if(all_obj_value[candidate[j]]<tmp_min)
            {
                tmp_min=all_obj_value[candidate[j]];
                candidate_best=candidate[j];
            }
        }
        
        new_pop.push_back(all_sol[candidate_best]);
    }
    
    return new_pop;
}
goaga::centroid goaga::vec_dim_transfer(population sol){
    centroid two_dim_sol=centroid(population_num,solution(cluster_num*dim));
    for(int i=0;i<population_num;i++)
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++)
                two_dim_sol[i][j*dim+k]=sol[i][j][k];
    return two_dim_sol;
}
goaga::population goaga::vec_dim_transfer(centroid sol){
    population three_dim_sol=population(population_num,centroid(cluster_num,solution(dim)));
    for(int i=0;i<population_num;i++)
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++)
                three_dim_sol[i][j][k]=sol[i][j*dim+k];
    return three_dim_sol;
}

//----------------------------------------------------------------------Crossover-------------------------------------------------------------------------------//
goaga::population goaga::crossover(centroid all_sol){
    int parent1,parent2,add_pos;
    centroid new_pop=centroid(population_num,solution(cluster_num*dim));
    while(all_sol.size()>0){
        
        if(all_sol.size()==1){
            new_pop[new_pop.size()-1]=all_sol[0];
            break;
        }
        parent1=rand()%all_sol.size();
        do{
            parent2=rand()%all_sol.size();
        }while(parent2==parent1);

        add_pos=population_num-all_sol.size();
        new_pop[add_pos]=all_sol[parent1];
        new_pop[add_pos+1]=all_sol[parent2];

        double r=(double)rand()/RAND_MAX;
        if(r<crossover_pro){
            int start_position=rand()%(dim*cluster_num-1);

            //做crossover
            for(int i=0;i<=start_position;i++){
                new_pop[add_pos][i]=all_sol[parent2][i];
                new_pop[add_pos+1][i]=all_sol[parent1][i];
            }
   
        }
        all_sol.erase(all_sol.begin()+parent1);
        all_sol.erase(all_sol.begin()+parent2);
    }
    return vec_dim_transfer(new_pop);
}
//----------------------------------------------------------------------Mutation--------------------------------------------------------------------------------//
void goaga::mutation(population& all_sol){
    for(int i=0;i<population_num;i++){
        double r=(double)rand()/RAND_MAX;
        if(r<mutation_pro){
            int mut_gene=rand()%(dim*cluster_num);
            r=((double)rand()/RAND_MAX)*2-1;
            all_sol[i][mut_gene/dim][mut_gene%dim]=all_sol[i][mut_gene/dim][mut_gene%dim]*r*2;
            bool flag=true;
            while(flag){
                flag=false;
                if(all_sol[i][mut_gene/dim][mut_gene%dim]<dmin)
                {
                    all_sol[i][mut_gene/dim][mut_gene%dim]+=(dmax-dmin)*((double)rand()/RAND_MAX);
                    flag=true;
                }
                else if(all_sol[i][mut_gene/dim][mut_gene%dim]>dmax)
                {
                   all_sol[i][mut_gene/dim][mut_gene%dim]-=(dmax-dmin)*((double)rand()/RAND_MAX);
                    flag=true;
                }
            } 
        }
    }
}
goaga::centroid goaga::run(){

    double tmp_best_obj_value;
    double avg_obj_value;
    population tmp_sol;
    vector<double>iter_obj_avg(numIter,0.0);
    
    init();
   
    for(int i=0;i<numRuns;i++){
        
        generate_new_solutions();
        
        fitness(currentSol);
        tmp_sol=currentSol;
        avg_obj_value=0.0;
        
        for(int j=0;j<numIter;j++){
            double choose_pro=0.01+0.98*((double)j/(double)numIter);
            double r=rand()/RAND_MAX;
            
            if((double)j>(double)numIter*0.75){
                double c=coef_c(j);
                for(int k=0;k<population_num;k++){
                    tmp_sol[k]=transition(currentSol,k,c);
                    currentSol[k]=tmp_sol[k];
                    /*
                    tmp_best_obj_value=evaluation(tmp_sol[k]);
                    
                    if(tmp_best_obj_value<best_obj_value)
                    {
                        best_obj_value=tmp_best_obj_value;
                        best_sol=tmp_sol[k];
                    }*/
                    
                }

                fitness(currentSol);
            }
            else
            {
                fitness(currentSol);
                currentSol=TS(currentSol);
                currentSol=crossover(vec_dim_transfer(currentSol));
                mutation(currentSol);
            }
            
            //cout<<"Iteration"<<j<<":"<<best_obj_value<<endl;
            iter_obj_avg[j]+=best_obj_value;
            
        }
         
    }
    
    ofstream plot_data("graph/plot_data/goaga_"+filename);
    for(int i=0;i<numIter;i++){
        cout<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;
        plot_data<<i<<" ";
        plot_data<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;
    }
    plot_data.close();

    return best_sol;

}


#endif
       