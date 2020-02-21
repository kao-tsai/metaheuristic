#ifndef __GOA2_H_INCLUDED__
#define __GOA2_H_INCLUDED__
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

class goa2{
	public:
		typedef vector<double> solution;
        typedef vector<solution> centroid;
        typedef vector<centroid> population;
		goa2(int,int,string,int,int,double,double);
		centroid run();
	private:
		void init();
        double dis_cal(solution,solution);
        void generate_new_solutions();
		double evaluation(centroid);
        double frand(double,double);
        double deg_forces(double);
        centroid transition(population,int,double);
        double coef_c(int);

        double normalize(double,int);

        void fitness(population&);

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
    centroid dim_minmax;
    double cmax;
    double cmin;
};

goa2::goa2(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
int x_cluster_num,
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
    cmax=xcmax;
    cmin=xcmin;
}

double goa2::deg_forces(double r){
    double f=0.5;
    double l=1.5;
    return f*exp((-r)/l)-exp(-r);
}

void goa2::fitness(population &sol) {
    double tmp_best;
	for (int i = 0; i < population_num; i++){
            tmp_best=evaluation(sol[i]);
            if(tmp_best<best_obj_value){
                    best_obj_value=tmp_best;
                    best_sol=sol[i];
            }
	}
}


void goa2::generate_new_solutions(){
    solution rand_centroid;
    solution tmp;
    for (int i = 0; i < population_num; i++)
        for(int j=0;j<cluster_num;j++)
            for(int k=0;k<dim;k++)
                currentSol[i][j][k]=frand(dim_minmax[k][0],dim_minmax[k][1]);
    
    best_obj_value=evaluation(currentSol[0]);
    best_sol=currentSol[0];
}
double goa2::evaluation(centroid sol){
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
double goa2::dis_cal(solution p1,solution p2){
    double sum=0.0;
    for(int i=0;i<dim;i++)
        sum+=pow(p1[i]-p2[i],2.0);
    return sum;
}

double goa2::coef_c(int cur_iter){
    double a;
    a=cmax-((cmax-cmin)/(double)(numIter-1))*(double)cur_iter;
    return a;
}


//----------------------------------------------------------------------force-------------------------------------------------------------------------------------//




goa2::centroid goa2::run(){

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
            currentSol=tmp_sol;
            
            fitness(currentSol);
            cout<<"Iteration"<<j<<":"<<best_obj_value<<endl;
            iter_obj_avg[j]+=best_obj_value;
            
        }
         
    }
    
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<(double)iter_obj_avg[i]/numRuns<<endl;
    
    
    return best_sol;

}


       
goa2::centroid goa2::transition(population sol,int sol_num,double cur_iter){
    double c,sum=0.0;
    double normalize_val;
    double tmp_dis;

    centroid tmp_sol=centroid(cluster_num,solution(dim));

    c=coef_c(cur_iter);
    
    for(int i=0;i<population_num;i++){
        if(sol_num==i)
            continue;
        
        for(int j=0;j<cluster_num;j++){
            tmp_dis=sqrt(dis_cal(sol[sol_num][j],sol[i][j]));
            if(tmp_dis==0)
                tmp_dis=0.00001;
            for(int k=0;k<dim;k++){
                normalize_val=normalize(fabs(sol[i][j][k]-sol[sol_num][j][k]),k);
                tmp_sol[j][k]+=c*((dim_minmax[k][1]-dim_minmax[k][0])/2.0)*deg_forces(normalize_val)*((sol[i][j][k]-sol[sol_num][j][k])/tmp_dis);
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
                if(tmp_sol[i][j]<dim_minmax[j][0])
                {
                   
                    tmp_sol[i][j]=tmp_sol[i][j]+(dim_minmax[j][1]-dim_minmax[j][0])*((double)rand()/RAND_MAX);
                    flag=true;
                }
                else if(tmp_sol[i][j]>dim_minmax[j][1])
                {
                    tmp_sol[i][j]= tmp_sol[i][j]-(dim_minmax[j][1]-dim_minmax[j][0])*((double)rand()/RAND_MAX);
                    flag=true;
                }
            }
            /*
            if(tmp_sol[i][j]<dim_minmax[j][0])
                tmp_sol[i][j]=dim_minmax[j][0];
            else if(tmp_sol[i][j]>dim_minmax[j][1])
                tmp_sol[i][j]=dim_minmax[j][1];
                */
        }

   return tmp_sol;
}

double goa2::normalize(double value,int wh_dim){
    return ((value-dim_minmax[wh_dim][0])/(dim_minmax[wh_dim][1]-dim_minmax[wh_dim][0]))*3+1;
}


double goa2::frand(double fmin,double fmax){
    double f = (double)(rand()%RAND_MAX) / RAND_MAX;
    return fmin + f * (fmax - fmin);
}


void goa2::init(){
    string line,con;
    solution tmp;
    int count_dim=0;
    if (!filename.empty()){
        ifstream file("GOA/"+filename);
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
            // getline(ss,con,',');
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
    currentSol=population(population_num,centroid(cluster_num,solution(dim)));
    
    }
    else{
		cout<<"Need Input!"<<endl;
    }
    
}

#endif