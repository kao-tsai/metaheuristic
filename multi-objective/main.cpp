#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include"evaluate.cpp"
#include"NSGAII/bnsgaii.h"
#include"SE/se.h"
// #include"MOPSO/mopso.h"
#include"MOPSO/mopso(noWH).h"
#include"SEDE/sede.h"
#include<fstream>
#include <math.h>
#include <algorithm>
using namespace std;

int main(int argc,char **argv)
{
    clock_t start,finish;
    double duration;
    // cout<<"Algorithm:"<<argv[1]<<endl;                      //nsga2
    // cout<<"Number of runs:"<<argv[2]<<endl;          //default 30
    // cout<<"Number of iterations:"<<argv[3]<<endl;//default 250
                 
   
    if(!strcmp(argv[1],"bnsga2"))
    {
        // cout<<"Number of population:"<<argv[4]<<endl;//default 500
        
        bnsgaii search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atof(argv[6]),argv[7]);
        
        search.run();

        int run_num=atoi(argv[2]);
        int pop_num=atoi(argv[4]);
        string func_name=argv[7];
        double tmp;
        double sum_IGD=0;
        double sum_SP=0;
        double sum_MS=0;
        ofstream output_igd;
        output_igd.open("pareto/"+func_name+"/bnsga2/IGD.txt");
        for(int i=0;i<run_num;i++){
            tmp=IGD(func_name,"bnsga2","opt_1000",i,2,pop_num,1000);
            sum_IGD+=tmp;
            output_igd<<tmp<<endl;
        }
        output_igd<<sum_IGD/run_num<<endl;
        output_igd.close();

        ofstream output_sp;
        output_sp.open("pareto/"+func_name+"/bnsga2/SP.txt");
        for(int i=0;i<run_num;i++){
            tmp=SP(func_name,"bnsga2","opt_1000",i,2,pop_num,1000);
            sum_SP+=tmp;
            output_sp<<tmp<<endl;
        }
        output_sp<<sum_SP/run_num<<endl;
        output_sp.close();
        
    }

    else if(!strcmp(argv[1],"mopso")){
        cout<<"Number of population:"<<argv[4]<<endl;     
        start=clock();
        mopso search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atoi(argv[6]),argv[7]);
        search.run();
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        
        int run_num=atoi(argv[2]);
        int pop_num=atoi(argv[6]);
        string func_name=argv[7];
        double tmp;
        double sum_IGD=0;
        double sum_SP=0;
        double sum_MS=0;
        ofstream output_igd;
        output_igd.open("pareto/"+func_name+"/mopso/IGD.txt");
        for(int i=0;i<run_num;i++){
            tmp=NSGAII_GD(func_name,"mopso","opt_1000",i,2,pop_num,1000);
            sum_IGD+=tmp;
            output_igd<<tmp<<endl;
        }
        output_igd<<sum_IGD/run_num<<endl;
        output_igd.close();

        ofstream output_sp;
        output_sp.open("pareto/"+func_name+"/mopso/SP.txt");
        for(int i=0;i<run_num;i++){
            tmp=NSGAII_SP(func_name,"mopso","opt_1000",i,2,pop_num,1000);
            sum_SP+=tmp;
            output_sp<<tmp<<endl;
        }
        output_sp<<sum_SP/run_num<<endl;
        output_sp.close();
        cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    }
    if(!strcmp(argv[1],"se"))
    {
        start=clock();
        se search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),argv[9],atof(argv[10]),atof(argv[11]));
        search.run();
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        int run_num=atoi(argv[2]);
        int pop_num=atoi(argv[4]);
        string func_name=argv[9];
        double tmp;
        double sum_IGD=0;
        double sum_SP=0;
        double sum_MS=0;
        ofstream output_igd;
        output_igd.open("pareto/"+func_name+"/se/IGD.txt");
        for(int i=0;i<run_num;i++){
            tmp=IGD(func_name,"se","opt_1000",i,2,pop_num,1000);
            sum_IGD+=tmp;
            output_igd<<tmp<<endl;
        }
        output_igd<<sum_IGD/run_num<<endl;
        output_igd.close();

        ofstream output_sp;
        output_sp.open("pareto/"+func_name+"/se/SP.txt");
        for(int i=0;i<run_num;i++){
            tmp=SP(func_name,"se","opt_1000",i,2,pop_num,1000);
            sum_SP+=tmp;
            output_sp<<tmp<<endl;
        }
        output_sp<<sum_SP/run_num<<endl;
        output_sp.close();
        cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    }
    if(!strcmp(argv[1],"sede"))
    {
        start=clock();
        sede search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),argv[9],atof(argv[10]),atof(argv[11]));
        search.run();
        finish=clock();
        
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        int run_num=atoi(argv[2]);
        int pop_num=atoi(argv[4]);
        string func_name=argv[9];
        double tmp;
        double sum_IGD=0;
        double sum_SP=0;
        double sum_MS=0;
        ofstream output_igd;
        output_igd.open("pareto/"+func_name+"/sede/IGD.txt");
        for(int i=0;i<run_num;i++){
            tmp=IGD(func_name,"sede","download_opt",i,2,pop_num,1000);
            sum_IGD+=tmp;
            output_igd<<tmp<<endl;
        }
        output_igd<<sum_IGD/run_num<<endl;
        output_igd.close();
        
        ofstream output_sp;
        output_sp.open("pareto/"+func_name+"/sede/SP.txt");
        for(int i=0;i<run_num;i++){
            tmp=SP(func_name,"sede","download_opt",i,2,pop_num,1000);
            sum_SP+=tmp;
            output_sp<<tmp<<endl;
        }
        output_sp<<sum_SP/run_num<<endl;
        output_sp.close();
        cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    }
    return 0;

}
