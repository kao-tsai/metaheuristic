#ifndef __ES_H_INCLUDED__
#define __ES_H_INCLUDE__
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<iomanip>
#include <time.h>
#include<vector>
#include<string>
#include<fstream>
using namespace std;

class hc{
        public:
                typedef  vector<int> solution;
                hc(int,int,int,string);
                solution run();

        private:
                void init();
                int evaluation(solution&);
                solution transition(solution);
                bool determination(int,int);
    private:
        int best_obj_value;
        int runs;
        int numPatterns;
        unsigned long long int iters;
        string filename;

        int f1 ;
        int f2;
        solution v;
        solution s;
        
};

hc::hc(int x_runs,int x_iters,int x_num_patterns,string x_filename)
{
        srand(time(0));
        runs=x_runs;
        iters=x_iters;
        numPatterns=x_num_patterns;
        filename=x_filename;
}
void hc::init(){
   
        if(!filename.empty()){
                ifstream file(filename);
                for(int i=0;i<numPatterns;i++)
                        file>>s[i];
                file.close();
        }
        else{
                s=solution(numPatterns);
               for(int i=0;i<numPatterns;i++)
                    s[i]=rand()%2;
        }
        
}

hc::solution hc::run()
{
        vector <double>avg_each_iter(iters,0);
        int avg_best_fitness=0;
      
         for(int i=0;i<runs;i++){
                init();
                f1=evaluation(s);
                for(int j=0;j<iters;j++)
                {      
                        v=transition(s);
                        f2=evaluation(v);
                        if(determination(f2,f1))
                        {
                                s=v;  
                                f1=f2;
                        }
                        avg_each_iter[j]+=f1;
                }
                avg_best_fitness+=f1;
        }
        for(int i=0;i<iters;i++)
                cout<<fixed<<setprecision(2)<<avg_each_iter[i]/runs<<endl;
        return  s;
}

hc::solution hc::transition(solution sol)
{       
        int i=rand()%numPatterns;
        sol[i]=!sol[i];
        return sol;
}

int hc::evaluation(solution& sol){
            int sum=0;

            for(int i=0;i<numPatterns;i++)
                   sum+=sol[i];

            return sum;
}
bool hc::determination(int a,int b)
{
        return a>b;                               
}
#endif