#ifndef __SA_H_INCLUDED__
#define __SA_H_INCLUDED__
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

class sa{
	public:
		typedef vector<int> solution;	
		sa(int,int,int,string,double,double);
		solution run();
	private:
		void init();
		int evaluation(solution&);
		solution transition(solution);
		bool determination(int,int,double);
		double annealing(double);

private:
	int numRuns;
	int numIter;
	int numPatterns;
	string filename;
	solution currentSol;
	solution tmpSol;
	int objValue;
	int tmpObjValue;
	int best_obj_value;
	solution bestSol;

	double minTemperature;
	double maxTemperature;
	double currentTemperature;
};
 
sa::sa(int xNumRuns,
int xNumIter,
int xNumPatterns,
string xfilename,
double xminTemperature,
double xmaxTemperature
)
{
    srand(time(0));
	numRuns=xNumRuns;
	cout<<numRuns<<endl;
	numIter=xNumIter;
	cout<<numIter<<endl;
	numPatterns=xNumPatterns;
	cout<<numPatterns<<endl;
	filename=xfilename;
	cout<<filename<<endl;
	minTemperature=xminTemperature;
	cout<<minTemperature<<endl;
	maxTemperature=xmaxTemperature;
	cout<<maxTemperature<<endl;
}

sa::solution sa::run(){
  double avg_obj_value;
  vector<double>iter_obj_avg(numIter,0.0);
  for(int i=0;i<numRuns;i++){
	init();
	avg_obj_value=0.0;
	objValue=evaluation(currentSol);
	for(int j=0;j<numIter;j++){
		//transition
		tmpSol=transition(currentSol);

		//evaluation
		tmpObjValue=evaluation(tmpSol);

		//determination

		if (objValue<tmpObjValue){
		    currentSol=tmpSol;
		    objValue=tmpObjValue;
		}
		else if(determination(tmpObjValue,objValue,currentTemperature))
		{
		    currentSol=tmpSol;
		    objValue=tmpObjValue;
		}
		if (best_obj_value<objValue)
		{
		    bestSol=currentSol;
		    best_obj_value=objValue;
		}
		currentTemperature=annealing(currentTemperature);
		iter_obj_avg[j]+=best_obj_value;
	}
	
  }

for(int i=0;i<numIter;i++)
	cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;

  return bestSol;
}
void sa::init(){
    currentSol=solution(numPatterns);
    objValue=0;
    best_obj_value=0;
    currentTemperature=maxTemperature;

    if (!filename.empty()){
        ifstream file(filename);
		for (int i=0;i<numPatterns;i++)
	    	file>>currentSol[i];
    }

    else{
		cout<<"randSol:";
		for (int i=0;i<numPatterns;i++)
		{
			currentSol[i]=rand()%2;
			
			cout<<currentSol[i];
		}
		cout<<endl;
    }
}

sa::solution sa::transition(solution sol){
	int i=rand()%numPatterns;
	sol[i]=!sol[i];
	return sol;
}
int sa::evaluation(solution& sol){
    int  count=0;
    for(int i=0;i<numPatterns;i++)
	count+=sol[i];
    return count;
}

bool sa::determination(int tmpObjValue,int objValue,double currentTemperature){
	double r=(double)rand()/(RAND_MAX);
	double p=exp((tmpObjValue-objValue)/currentTemperature);
	return r<p;
}

double sa::annealing(double temperature){
	return 0.95*temperature;
}

#endif
