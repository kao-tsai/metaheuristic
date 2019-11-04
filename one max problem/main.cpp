#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"GA/ga.h"
#include"ACO/aco.h"
using namespace std;


int main(int argc,char **argv)
{
    cout<<"Algorithm:"<<argv[1]<<endl;
    cout<<"Number of runs:"<<argv[2]<<endl;
    cout<<"Number of iterations:"<<argv[3]<<endl;
    cout<<"Name of file:"<<argv[4]<<endl;
    cout<<"Number of population:"<<argv[5]<<endl;
    cout<<"alpha="<<argv[6]<<endl;
    cout<<"beta="<<argv[7]<<endl;
    cout<<"rho="<<argv[8]<<endl;
    
    if(!strcmp(argv[1],"ga"))
    {
        ga search(argv[2],atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6],atoi(argv[7]));
        ga::solution sol=search.run();
        
        for(int i=0;i<sol.size();i++)
            cout<<sol[i];
        cout<<endl;
    }

    if(!strcmp(argv[1],"aco"))
    {
        //run  iter  filename  ants  alpha  beta  rho  mod
        aco search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8]),atoi(argv[9]));
        aco::solution sol=search.run();
        for(int i=0;i<sol.size();i++){
             if(i==sol.size()-1){
                cout<<endl<<"best_obj_value:"<<fixed<<setprecision(3)<<sol[i]<<endl;    
                break;
             }
            cout<<fixed<<setprecision(0)<<sol[i]<<" ";
           
        }
        
        

    }
    return 0;
}