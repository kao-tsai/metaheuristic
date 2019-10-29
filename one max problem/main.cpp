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
    cout<<"Number of patterns:"<<argv[4]<<endl;
    cout<<"Name of file:"<<argv[5]<<endl;
    cout<<"Number of population:"<<argv[6]<<endl;
    
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
       
        aco search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5],atoi(argv[6]));
        aco::solution sol=search.run();
        //cout<<sol.size()<<"ehy";
        for(int i=0;i<sol.size();i++){
            cout<<sol[i];
            if(i==sol.size()-1);
                break;
            cout<<" ";
        }
        cout<<endl;

    }
    return 0;
}