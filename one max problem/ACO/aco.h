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
		aco(int,int,int,string,int);
		solution run();
	private:
		void init();
		int evaluation(solution&);
        void ant_activity();
        double dis_cal(double,double,double,double);
		void phe_update();
private:
	int numRuns;
	int numIter;
	int numPatterns;
	string filename;

	population currentSol;
	solution obj_value;
	int best_obj_value;
	solution best_sol;
    int ant_num;
    int city_num;
    population node;
	population node_dis;
    population pheromone;
	population add_phe;
	double rho;//揮發率
    int alpha,beta;
};
aco::aco(int xNumRuns,
int xNumIter,
int xNumPatterns,
string xfilename,
int xPopulationNum
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
	numPatterns=xNumPatterns;
	filename=xfilename;
    ant_num=xPopulationNum;
	rho=0.9;
}

double aco::dis_cal(double x1,double y1,double x2,double y2){
	
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));

}



void aco::init(){
	best_obj_value=1000;
    //objvalue_reset
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
	
    city_num=node.size();
	//cout<<city_num<<"this";
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
                    pheromone[i][j]=0.00017;
                    pheromone[j][i]=0.00017;
                }

                else{
                    node_dis[i][j]=0;
                    pheromone[i][j]=0;
					
                }            
        }
    alpha=1;
    beta=2;
    }
    else{
		cout<<"Need Input!"<<endl;
    }
	
}



void aco::ant_activity(){
    solution ants(city_num);
	solution prob(city_num-1);
	int set,stepOn;
	double sum=0,r,pathLength=0;
	
	for (int i = 0; i < city_num; i++)
			add_phe[i].resize(city_num,0.0);
	
	
	for (int i = 0; i < ant_num; i++)
	{
		set = city_num-1;
		ants.resize(city_num);
		stepOn = i;
		for (int j = 0; j < city_num; j++)
			ants[j]=j;

		ants.erase(ants.begin()+stepOn);
		ants.push_back(stepOn);

		pathLength = 0;
		for (int j = 0; j < city_num; j++)
		{
			sum = 0;
			for (int k = 0; k < set;k++)
			{
				prob[k]=pow(pheromone[stepOn][ants[k]],alpha) * pow((1.0 / node_dis[stepOn][ants[k]]), beta);
				sum+=prob[k];

			}
			cout<<sum<<"this";
			for (int k = 0; k < set; k++)
				if (k == 0)
					prob[k] = prob[k] / sum;
				else
					prob[k] = prob[k - 1] + (prob[k] / sum);
			
			r = (double)rand() / (RAND_MAX);
			for (int m = 0; m < set;m++)
				if (r < prob[m])
				{
					cout<<ants[m]<<endl;
					pathLength += node_dis[stepOn][ants[m]];
					stepOn = ants[m];
					ants.push_back(ants[m]);
					ants.erase(ants.begin() + m);
					break;
				}
			set--;
			
		}

		ants.push_back(ants[0]);
		pathLength += node_dis[ants[city_num - 1]][ants[0]];
		ants.push_back(pathLength);
		for (int i = 0; i < city_num; i++)
			add_phe[i].resize(city_num,0.0);

		for (int j = 0; j < city_num; j++)
			if (ants[j] > ants[j + 1])
			{
				add_phe[ants[j + 1]][ants[j]] += 100 / pathLength;
				add_phe[ants[j]][ants[j + 1]] = add_phe[ants[j + 1]][ants[j]];
			}
			else
			{
				add_phe[ants[j]][ants[j + 1]] += 100 / pathLength;
				add_phe[ants[j + 1]][ants[j]] = add_phe[ants[j]][ants[j + 1]];
			}

		if (pathLength < best_obj_value)
		{
			best_obj_value = pathLength;
			best_sol = ants;
		}
	
	}
}


void aco::phe_update() {
	double temp;
	for (int i = 0; i < city_num; i++)
		for (int j = i; j < city_num; j++)
			if (i == j)
				continue;
			else{
				pheromone[i][j] = (1 - rho) * pheromone[i][j] + add_phe[i][j];
				pheromone[j][i] =pheromone[i][j];
			}
}

aco::solution aco::run(){
	
    double avg_obj_value;
    vector<double>iter_obj_avg(numIter,0.0);
    for(int i=0;i<numRuns;i++){
        init();
        avg_obj_value=0.0;
		cout<<i<<endl;
        for(int j=0;j<numIter;j++){
            ant_activity();
			phe_update();
			
            iter_obj_avg[j]+=best_obj_value;
        }
        
    }

    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;

    return best_sol;

/*
for(int i=0;i<city_num;i++)
{
	for(int j=0;j<city_num;j++)
		cout<<node_dis[i][j]<<" ";
	cout<<endl;
}
return best_sol;*/
}
#endif