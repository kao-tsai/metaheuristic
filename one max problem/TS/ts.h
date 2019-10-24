#ifndef __TS_H_INCLUDED__
#define __TS_H_INCLUDED__
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

class ts{
	public:
		typedef vector<int> solution;	
		ts(int,int,int,string,int);
		solution run();
	private:
		void init();
		int evaluation(solution&);
		solution transition(solution);
		bool determination(int,int);
		solution noInTabuNeighbor(solution);
		bool checkInTabu(solution&);
		void push_in_tabu(solution sol);
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
	vector<string> tabuList;
	int tabu_list_max_size;
};
 
ts::ts(int xNumRuns,
int xNumIter,
int xNumPatterns,
string xfilename,
int xtabuSize
)
{
    srand(time(0));
	numRuns=xNumRuns;
	numIter=xNumIter;
	numPatterns=xNumPatterns;
	filename=xfilename;
	tabu_list_max_size=xtabuSize;
}

ts::solution ts::run(){
  double avg_obj_value;
  vector<double>iter_obj_avg(numIter,0.0);
  for(int i=0;i<numRuns;i++){
	init();
	objValue=evaluation(currentSol);
	for(int j=0;j<numIter;j++){
			//transition
			tmpSol=noInTabuNeighbor(currentSol);

			//evaluation
			tmpObjValue=evaluation(tmpSol);

			//determination
			if(determination(tmpObjValue,objValue))
			{
					currentSol=tmpSol;
					objValue=tmpObjValue;
					if(objValue>best_obj_value)
					{
							best_obj_value=objValue;
							bestSol=currentSol;
					}
					push_in_tabu(currentSol);
			}
			
			iter_obj_avg[j]+=best_obj_value;
	}
	
  }

for(int i=0;i<numIter;i++)
	cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;

  return bestSol;
}
void ts::init(){
    currentSol=solution(numPatterns);
    objValue=0;
    best_obj_value=0;
	tabuList.resize(0);
    if (!filename.empty()){
        ifstream file(filename);
		for (int i=0;i<numPatterns;i++)
	    	file>>currentSol[i];
    }

    else{
		for (int i=0;i<numPatterns;i++)
				currentSol[i]=rand()%2;

    }
}

ts::solution ts::noInTabuNeighbor(solution sol){
		solution tmp;
		int count=0;
		tmp=transition(sol);
		while(1)
		{
			if(!checkInTabu(tmp))			
							break;
			
			tmp=transition(sol);
		}
		
		return tmp;

}
bool ts::checkInTabu(solution& sol){
	string buffer="";
	stringstream ss;
	for(int i=0;i<sol.size();i++)
		 ss<<sol[i];
	buffer=buffer+ss.str();
	for(int i=0;i<tabuList.size();i++)
			if(tabuList[i]==buffer)
					return true;

	return false;
}
void ts::push_in_tabu(solution sol){
		string buffer="";
		stringstream ss;
		for(int i=0;i<sol.size();i++)
				 ss<<sol[i];
				
		buffer=buffer+ss.str();
		tabuList.push_back(buffer);
		if(tabuList.size()>tabu_list_max_size)
				tabuList.erase(tabuList.begin());
	
}
ts::solution ts::transition(solution sol){
	int i=rand()%numPatterns;
	sol[i]=!sol[i];
	return sol;
}
int ts::evaluation(solution& sol){
    int  count=0;
    for(int i=0;i<numPatterns;i++)
	count+=sol[i];
    return count;
}

inline bool ts::determination(int f2,int f1){
	return f2>f1;
}


#endif
