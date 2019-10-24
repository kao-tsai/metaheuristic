#ifndef __SA_H_INCLUDED__
#define __SA_H_INCLUDED__
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<iomanip>
#include <time.h>
#include<vector>
#include<string>
#include<fstream>
#include<math.h>
using namespace std;

class sa{
        public:
                typedef  vector<int> solution;
                sa(int,int,int,string,double,double);
                solution run();

        private:
                void init();
                int evaluation(solution&);
                solution transition(solution);
                bool determination(int,int,double);
                bool checkCapacity(solution&);
                double annealing(double);
    private:
        int best_obj_value;
        int runs;
        int numPatterns;
        unsigned long long int iters;
        string filename;

        int objValue ;
        int tmpObjValue;
        solution bestSol;
        solution tmpSol;
        solution weight;
        solution value;
        solution currentSol;
        int capacity;
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
	runs=xNumRuns;
	iters=xNumIter;
	numPatterns=xNumPatterns;
	filename=xfilename;
	minTemperature=xminTemperature;
	maxTemperature=xmaxTemperature;
        weight=solution(numPatterns);
        value=solution(numPatterns);

}
void sa::init(){
        objValue=0;
        best_obj_value=0;
        currentTemperature=maxTemperature;
        if(!filename.empty()){

                ifstream file(filename);
                file>>capacity;
                for(int i=0;i<numPatterns;i++)     
                        file>>weight[i];
                for(int i=0;i<numPatterns;i++)     
                        file>>value[i];
                file.close();

                currentSol=solution(numPatterns);
                for(int i=0;i<numPatterns;i++)
                        currentSol[i]=rand()%2;
                while(checkCapacity(currentSol))
                {
                        currentSol.resize(numPatterns,0);
                        for(int i=0;i<numPatterns;i++)
                                currentSol[i]=rand()%2;
                }
        }
        else{
                currentSol=solution(numPatterns);
               for(int i=0;i<numPatterns;i++)
                    currentSol[i]=rand();
        }	
}

sa::solution sa::run()
{
        vector <double>avg_each_iter(iters,0);
        int avg_best_fitness=0;
        int mar=100;
         for(int i=0;i<runs;i++){
                init();
                objValue=evaluation(currentSol);
                for(int j=0;j<iters;j++)
                {      
                        for(int k=0;k<mar;k++){
                       	//transition
		        tmpSol=transition(currentSol);

		        //evaluation
		        tmpObjValue=evaluation(tmpSol);
                      
                        if (objValue<tmpObjValue){
                                currentSol=tmpSol;
                                objValue=tmpObjValue;
                        }
                        else if(determination(tmpObjValue,objValue,currentTemperature)){
                                currentSol=tmpSol;
                                objValue=tmpObjValue;
                        }
                        if (best_obj_value<objValue){
                                bestSol=currentSol;
                                best_obj_value=objValue;
                        }
                        }
                        currentTemperature=annealing(currentTemperature);
                        avg_each_iter[j]+=best_obj_value;
                }
                avg_best_fitness+=best_obj_value;
        }
        for(int i=0;i<iters;i++)
                cout<<fixed<<setprecision(2)<<avg_each_iter[i]/runs<<endl;
        return  bestSol;
}
bool sa::checkCapacity(solution& sol){
        int count=0;
        for(int i=0;i<numPatterns;i++)
        {
                if(sol[i]==1)
                        count+=weight[i];
                if(count>capacity)
                        return true;
        }
        return false;
}
sa::solution sa::transition(solution sol)
{         /*
        int first,second;

	first = rand() % numPatterns;	
	sol[first] = !sol[first];

	while (1) {
		if (checkCapacity(sol))
		{
			do {
				second = rand() % numPatterns;
			} while (sol[second] == 0 && second == first);

			sol[second] = 0;
		}
		else
			break;
	}
        */
         
        int i=rand()%numPatterns;
        sol[i]=!sol[i];
        while(checkCapacity(sol))
        {
               sol[i]=!sol[i];
               i=rand()%numPatterns;
               sol[i]=!sol[i];
        }
        /*
        bool first;
        first=!checkCapacity(sol);
        while(checkCapacity(sol) || first)
        {
                i=rand()%numPatterns;
                while(!sol[i])
                        i=rand()%numPatterns;
                sol[i]=!sol[i];

                i=rand()%numPatterns;
                while(sol[i])
                        i=rand()%numPatterns;
                sol[i]=!sol[i];
                first=false;
        }*/
        return sol;
}

int sa::evaluation(solution& sol){
            int sum=0;

            for(int i=0;i<numPatterns;i++)
                   if(sol[i])
                        sum+=value[i];

            return sum;
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