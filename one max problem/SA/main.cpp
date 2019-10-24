#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"sa.h"
using namespace std;


int main(int argc,char **argv)
{
    cout<<"Algorithm:"<<argv[1]<<endl;
    cout<<"Number of runs:"<<argv[2]<<endl;
    cout<<"Number of iterations:"<<argv[3]<<endl;
    cout<<"Number of patterns:"<<argv[4]<<endl;
    cout<<"Name of file:"<<argv[5]<<endl;
    if(!strcmp(argv[1],"sa"))
    {
        sa search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5],atof(argv[6]),atof(argv[7]));
        sa::solution sol=search.run();
        
        for(int i=0;i<sol.size();i++)
            cout<<sol[i];
        cout<<endl;
    }
    
    return 0;
}

