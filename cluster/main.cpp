#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"GOA/goa.h"
#include"GOA/goa2.h"
using namespace std;


int main(int argc,char **argv)
{
    cout<<"Algorithm:"<<argv[1]<<endl;
    cout<<"Number of runs:"<<argv[2]<<endl;
    cout<<"Number of iterations:"<<argv[3]<<endl;
    cout<<"Name of file:"<<argv[4]<<endl;
    cout<<"Number of population:"<<argv[5]<<endl;
    /*
    if(!strcmp(argv[1],"ga"))
    {
        cout<<"Crossover Rate:"<<argv[6]<<endl;
        cout<<"Mutation Rate:"<<argv[7]<<endl;
        cout<<"Selection:"<<argv[8]<<endl;
        cout<<"Crossover Method:"<<argv[9]<<endl;

        ga search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atof(argv[6]),atof(argv[7]),argv[8],argv[9]);
        ga::solution sol=search.run();
        for(int i=0;i<sol.size();i++)
            cout<<sol[i]<<" ";
        cout<<endl;
    }*/
    if(!strcmp(argv[1],"goa"))
    {
        cout<<"Number of clusters:"<<argv[6]<<endl;
        cout<<"cmax:"<<argv[7]<<endl;
        cout<<"cmin:"<<argv[8]<<endl;
     

        goa search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atoi(argv[6]),atof(argv[7]),atof(argv[8]));
        goa::centroid sol=search.run();
        
        
        for(int i=0;i<sol.size();i++){
            cout<<"Cluster"<<i<<":";
            for(int j=0;j<sol[i].size();j++)
                cout<<sol[i][j]<<" ";
            cout<<endl;
        }
       
    }

    if(!strcmp(argv[1],"goa2"))
    {
        cout<<"Number of clusters:"<<argv[6]<<endl;
        cout<<"cmax:"<<argv[7]<<endl;
        cout<<"cmin:"<<argv[8]<<endl;
     

        goa2 search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atoi(argv[6]),atof(argv[7]),atof(argv[8]));
        goa2::centroid sol=search.run();
        
        
        for(int i=0;i<sol.size();i++){
            cout<<"Cluster"<<i<<":";
            for(int j=0;j<sol[i].size();j++)
                cout<<sol[i][j]<<" ";
            cout<<endl;
        }
       
    }
    return 0;
}