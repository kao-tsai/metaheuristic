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
    
    if(!strcmp(argv[1],"goa"))
    {
        cout<<"cmax:"<<argv[6]<<endl;
        cout<<"cmin:"<<argv[7]<<endl;
        cout<<"Dimension:"<<argv[8]<<endl;

        goa search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atof(argv[6]),atof(argv[7]),atoi(argv[8]));
        goa::solution sol=search.run();
        
        


        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
        
    }

    if(!strcmp(argv[1],"goa2"))
    {
        cout<<"cmax:"<<argv[6]<<endl;
        cout<<"cmin:"<<argv[7]<<endl;
        cout<<"Dimension:"<<argv[8]<<endl;

        goa search(atoi(argv[2]),atoi(argv[3]),argv[4],atoi(argv[5]),atof(argv[6]),atof(argv[7]),atoi(argv[8]));
        goa::solution sol=search.run();

        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
        
       
    }

    return 0;
}