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
		ga(string,int,int,int,string,int);
		solution run();
	private:
		void init();
		int evaluation(solution&);
		population crossover(population);
        population mutation(population);
        solution fitness(population);
		population RWS(population);
        population TS(population);

private:
	int numRuns;
	int numIter;
	int numPatterns;
	string filename;
    string selection;
	population currentSol;
	solution obj_value;
	int best_obj_value;
	solution best_sol;
    int population_num;
	double mutation_pro;
    double crossover_pro;
};
ga::ga(string xSelection,int xNumRuns,
int xNumIter,
int xNumPatterns,
string xfilename,
int xPopulationNum
)
{
    srand(time(0));
    selection=xSelection;
    numRuns=xNumRuns;
	numIter=xNumIter;
	numPatterns=xNumPatterns;
	filename=xfilename;
    crossover_pro=0.6;
    mutation_pro=0.1;
    population_num=xPopulationNum;
}

ga::solution ga::fitness(population all_sol) {
    solution tmp(population_num);
	for (int i = 0; i < population_num; i++){
            tmp[i]=evaluation(all_sol[i]);

            if(best_obj_value<tmp[i]){
                    best_obj_value=tmp[i];
                    best_sol=all_sol[i];
            }
	}
    return tmp;
}
void ga::init(){
    best_obj_value=0;
    currentSol.clear();
    currentSol.resize(population_num);
	for (int i = 0; i < population_num; i++)
        for(int j=0;j<numPatterns;j++)
            currentSol[i].push_back(rand()%2);

}
ga::population ga::RWS(population all_sol){
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


ga::population ga::TS(population all_sol){
    population new_pop;
    bool flag;
    int tournament_size=3,tmp,max,candidate_best;
    solution candidate;
    for(int i=0;i<population_num;i++){

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
        

        max=evaluation(all_sol[candidate[0]]);
        candidate_best=candidate[0];
        for(int j=1;j<tournament_size;j++){
            tmp=evaluation(all_sol[candidate[j]]);
            if(tmp>max)
            {
                max=tmp;
                candidate_best=candidate[j];
            }

        }
        new_pop.push_back(all_sol[candidate_best]);
    }
    return new_pop;
}

ga::population ga::crossover(population all_sol){
        int cut,parent1,parent2;
        population new_pop;
        solution temp_p1;
        solution temp_p2;
         solution temp(numPatterns);

        for(int i=0;i<population_num;i=i+2){

            cut=rand()%(numPatterns-1);
            temp_p1=all_sol[i];
            temp_p2=all_sol[i+1];
            temp=temp_p1;

            //crossover probability
            double r=(double)rand()/RAND_MAX;
            if(r<crossover_pro)
                for (int j = 0; j < cut; j++){
                    temp_p1[j] = temp_p2[j];
                    temp_p2[j] = temp[j];
                }
            
            if(i==population_num-1)
                new_pop.push_back(all_sol[population_num-1]);
            
            else{
                new_pop.push_back(temp_p1);
                new_pop.push_back(temp_p2);
            }
        }
        return new_pop;
}

ga::population ga::mutation(population all_sol){
    int mutation_position;
    for(int i=0;i<population_num;i++){
            double r=(double)rand()/RAND_MAX;
            if(r<mutation_pro){
                mutation_position=rand()%numPatterns;
                all_sol[i][mutation_position]=!all_sol[i][mutation_position];
            }
    }
    return all_sol;
}

ga::solution ga::run(){
    int count=0;
    double avg_obj_value;
    vector<double>iter_obj_avg(numIter,0.0);
    for(int i=0;i<numRuns;i++){
        //cout<<"hello:"<<count<<endl;
        count++;
        init();
        avg_obj_value=0.0;
        for(int j=0;j<numIter;j++){
            
            obj_value=fitness(currentSol);
            if(selection=="ts")
                currentSol=TS(currentSol);
            else if(selection=="rws")
                currentSol=RWS(currentSol);
            
            currentSol=crossover(currentSol);
            currentSol=mutation(currentSol);

            iter_obj_avg[j]+=best_obj_value;
        }
        
    }

    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;

    return best_sol;

}

int ga::evaluation(solution& sol){
    int  count=0;
    for(int i=0;i<numPatterns;i++)
	count+=sol[i];
    return count;
}
#endif
       