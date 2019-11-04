#ifndef __ACO_H_INCLUDED__
#define __ACO_H_INCLUDED__
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include<iomanip>
using namespace std;


class aco{
	public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		aco(int,int,string,int,double,double,double,int);
		solution run();
	private:
		void init();
		int evaluation(solution&);
        void ant_activity();
		void reset_phe_table();
		void two_opt();
		double greedy_search();
        double dis_cal(double,double,double,double);
		void phe_update();
private:
	int numRuns;
	int numIter;
	int numPatterns;
	string filename;

	double all_run_best_obj;
	solution best_sol;
    int ant_num;
    int city_num;
	double init_phe;
    population node;
	population node_dis;
    population pheromone;
	population add_phe;
	double each_run_best_obj;
	int mod;
	double rho;//揮發率
	double alpha,beta;
};
aco::aco(int xNumRuns,
int xNumIter,
string xfilename,
int xPopulationNum,
double xalpha,
double xbeta,
double xrho,
int xmod
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	filename=xfilename;
    ant_num=xPopulationNum;
	alpha=xalpha;
	beta=xbeta;
	rho=xrho;
	mod=xmod;
}

double aco::dis_cal(double x1,double y1,double x2,double y2){
	
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));

}



void aco::init(){
    string line,con;
    solution tmp;
    if (!filename.empty()){
        ifstream file(filename);
		node.clear();
		while(file>>line){
            file>>line;
            tmp.push_back(atof(line.c_str()));
            file>>line;
            tmp.push_back(atof(line.c_str()));
            node.push_back(tmp);
            tmp.clear();
        }
	file.close();
    city_num=node.size();

    node_dis=population(city_num);
    pheromone=population(city_num);
	add_phe=population(city_num);
	for(int i=0;i<city_num;i++){
			node_dis[i].resize(city_num);
			pheromone[i].resize(city_num);
			
	}
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
                    pheromone[i][j]=0;
				
                }            
        }
	init_phe=greedy_search();
	for(int i=0;i<city_num;i++)
        for(int j=i;j<city_num;j++)
                if(i!=j){
					
                    pheromone[i][j]=init_phe;
                    pheromone[j][i]=init_phe;
                }

				
                     
    
    }
    else{
		cout<<"Need Input!"<<endl;
    }
	
}

void aco::reset_phe_table(){
	  for(int i=0;i<city_num;i++)
        for(int j=i;j<city_num;j++)
        {		
                if(i!=j){
                    pheromone[i][j]=init_phe;
                    pheromone[j][i]=init_phe;	
                }

                else
                    pheromone[i][j]=0;				
            
        }
}

double aco::greedy_search(){
	int start_city,min_city,first_city;
	double pathLength=0.0;
	double min_dis;
	vector<int>remain_city(city_num);
	//best_sol.clear();
	start_city=rand()%city_num;
	first_city=start_city;
	for(int i=0;i<city_num;i++)
		remain_city[i]=i;

	//best_sol.push_back(remain_city[start_city]);

	remain_city.erase(remain_city.begin()+start_city);

	while(1)
	{
		if(remain_city.size()==0)
		{
			pathLength+=node_dis[first_city][start_city];
			//best_sol.push_back(first_city);
			//best_sol.push_back(pathLength);
			
			return 1.0/(pathLength*city_num);
		}
		min_dis=node_dis[remain_city[0]][start_city];
		min_city=0;
		for(int j=0;j<remain_city.size();j++)
			if(min_dis<node_dis[start_city][remain_city[j]]){	
				min_city=j;
				min_dis=node_dis[start_city][remain_city[j]];
			}
		
		start_city=remain_city[min_city];
		//best_sol.push_back(start_city);
		pathLength+=min_dis;

		remain_city.erase(remain_city.begin()+min_city);
	}
	
}

void aco::ant_activity(){
    solution remain_city(city_num);
	solution prob(city_num-1);
	solution reset_remain(city_num);
	population ants(ant_num);
	int stepOn;
	double sum=0,r,pathLength=0;
	for(int i=0;i<city_num;i++)
		reset_remain[i]=i;
	for (int i = 0; i < city_num; i++){
		add_phe[i].clear();
		add_phe[i].resize(city_num,0.0);
	
	}

	for (int i = 0; i < ant_num; i++)
	{
		//start city
		stepOn = rand()%city_num;

		remain_city=reset_remain;
		remain_city.erase(remain_city.begin()+stepOn);
		ants[i].push_back(stepOn);

		pathLength = 0;
		
		while(remain_city.size()>0)
		{
			sum = 0;

			for (int k = 0; k < remain_city.size();k++)
			{
				
				prob[k]=pow(pheromone[stepOn][remain_city[k]],alpha) * pow((1.0 / node_dis[stepOn][remain_city[k]]),beta);
				sum+=prob[k];
				//cout<<k<<endl;
				
			}
	
			for (int k = 0; k < remain_city.size(); k++)
				if (k == 0){
					prob[k] = prob[k] / sum;
					//cout<<k<<":"<<prob[k]<<endl;
				}
				else{
					prob[k] = prob[k - 1] + (prob[k] / sum);
					//cout<<k<<":"<<prob[k]<<endl;
				}
					
				

			r = (double)rand() / (RAND_MAX);
			
			for (int m = 0; m < remain_city.size();m++)
				if (r < prob[m])
				{
					// if(m==20)
					// 	cout<<"hello"<<endl;
					pathLength += node_dis[stepOn][remain_city[m]];
					stepOn = remain_city[m];
					ants[i].push_back(remain_city[m]);
					remain_city.erase(remain_city.begin() + m);
					break;
				}
			
			
		}
		
		ants[i].push_back(ants[i][0]);
		pathLength += node_dis[ants[i][city_num-1]][ants[i][0]];
		ants[i].push_back(pathLength);
		
		for (int j = 0; j < city_num-1; j++){
			//if (ants[i][j] > ants[i][j + 1])
				add_phe[ants[i][j + 1]][ants[i][j]] +=1/ pathLength;
				add_phe[ants[i][j]][ants[i][j+1]] =add_phe[ants[i][j + 1]][ants[i][j]];
			//else
			//	add_phe[ants[i][j]][ants[i][j + 1]] += 1/ pathLength;	
			
		}
		
		if (pathLength < each_run_best_obj)
		{
			each_run_best_obj = pathLength;
			best_sol = ants[i];
		}
		
	}
	
}

void aco::two_opt(){
	int reverse1,reverse2,tmp;
	solution init_best_sol;
	solution opt_best_sol;
	init_best_sol=best_sol;
	opt_best_sol=best_sol;
	
	double pathLength=0.0,opt_best_obj;
	opt_best_obj=each_run_best_obj;
	for(int i=1;i<city_num-1;i++){
		for(int j=i+1;j<city_num;j++)
		{
			
			pathLength=0;
			reverse1=i;
			reverse2=j;
			
			while(reverse1<reverse2)
			{
				tmp=best_sol[reverse1];
				best_sol[reverse1]=best_sol[reverse2];
				best_sol[reverse2]=tmp;
				reverse1++;
				reverse2--;
			}
			
			for(int i=0;i<best_sol.size()-2;i++)
				pathLength+=node_dis[best_sol[i]][best_sol[i+1]];
			
			if(pathLength<opt_best_obj){
				opt_best_obj=pathLength;
				opt_best_sol=best_sol;
				opt_best_sol.erase(opt_best_sol.end()-1);
				
				opt_best_sol.push_back(pathLength);
			}
			best_sol=init_best_sol;
			
		}
	}
	best_sol=opt_best_sol;
	each_run_best_obj=opt_best_obj;
	
}
void aco::phe_update() {
	double temp;
	for (int i = 0; i < city_num; i++)
		for (int j = i; j < city_num; j++)
			if (i != j){
				//if(add_phe[i][j]!=0)
				pheromone[i][j] = (1 - rho) *pheromone[i][j]+add_phe[i][j];
				if(pheromone[i][j]<init_phe)
				 	pheromone[i][j]=init_phe;
				 
				pheromone[j][i] =pheromone[i][j];
			}
	
}

aco::solution aco::run(){
	
    double avg_obj_value;
	clock_t start,finish;
    solution iter_obj_avg(numIter,0.0);
	solution all_run_best_sol;
	init();
	all_run_best_obj=2000;
	start=clock();
    for(int i=0;i<numRuns;i++){
        reset_phe_table();
		each_run_best_obj=2000;
        avg_obj_value=0.0;

        for(int j=0;j<numIter;j++){
			
            ant_activity();
			
			phe_update();
			
			two_opt();
			// if(j==numIter-1)
			// {
			// 	for (int i = 0; i < city_num; i++){
			// 			for (int j = 0; j < city_num; j++)
			// 				cout<<pheromone[i][j]<<" ";
			// 			cout<<endl;
			// 	}
			// }
            iter_obj_avg[j]+=each_run_best_obj;
			
        }
        if(each_run_best_obj<all_run_best_obj)
		{
			all_run_best_obj=each_run_best_obj;
			all_run_best_sol=best_sol;
		}
		

    }
	finish=clock();
	if(mod==0){
		ofstream output1("ACO/alpha_"+to_string(alpha)+"_plot_data.txt");
		for(int i=0;i<numIter;i++){
			cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;	
			output1<<i<<" "<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;	
		}
		output1.close();
	}
	else if(mod==1){
		ofstream output2("ACO/beta_"+to_string(beta)+"_plot_data.txt");
		for(int i=0;i<numIter;i++){
			cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
			output2<<i<<" "<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;			

		}
		output2.close();
	}
	else if(mod==3){
		ofstream output3("ACO/rho_"+to_string(rho)+"_plot_data.txt");
		for(int i=0;i<numIter;i++){
			cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
			output3<<i<<" "<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;			

		}
		output3.close();
	}
	else{
		ofstream output4("ACO/best_data.txt");
		for(int i=0;i<numIter;i++){
			cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;
			output4<<i<<" "<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;			

		}
		output4.close();

	}
    
	
	cout << "wholeTime:" << (double)(finish- start) / CLOCKS_PER_SEC<<endl;

	ofstream cout_path("ACO/path_data.txt");

	for(int i=0;i<all_run_best_sol.size()-1;i++)
			cout_path<<node[all_run_best_sol[i]][0]<<" "<<node[all_run_best_sol[i]][1]<<endl;

	
	
	

    return all_run_best_sol;
	
}
#endif