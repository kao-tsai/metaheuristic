#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"HHO/HHO.h"
#include"SE/SE.h"
#include"MFO/mfo.h"
#include"PSO/pso.h"
#include"MFODE/mfode.h"
#include"GA/ga.h"
#include"TLBO/tlbo.h"
#include<algorithm>
#include<functional>
using namespace std;


int main(int argc,char **argv)
{
    // cout<<"Algorithm:"<<argv[1]<<endl;
    // cout<<"Number of runs:"<<argv[2]<<endl;
    // cout<<"Number of iterations:"<<argv[3]<<endl;
    
    
    if(!strcmp(argv[1],"hoo"))
    {

        cout<<"Number of population:"<<argv[4]<<endl;
        hoo search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
        hoo::solution sol=search.run();
        
        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
        /*
        clock_t start,finish;
        double duration;
        vector<double> a(1000000,-1.2985854437);
        vector<double> b(1000000,245.4522);
        vector<double> c(1000000,-4.5);
        // transform(a.begin(),a.end(),c.begin(),static_cast<double(*)(double)>(&fabs));
        start=clock();
        transform(a.begin(),a.end(),b.begin(),c.begin(),multiplies<double>());
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        cout<<"transform time:"<<duration<<endl;
        start=clock();
        for(int i=0;i<1000000;i++)
           c[i]=a[i]+b[i];
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        cout<<"Loop time:"<<duration<<endl;
        */
    }
    if(!strcmp(argv[1],"tlbo"))
    {

        cout<<"Number of population:"<<argv[4]<<endl;
        tlbo search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
        tlbo::solution sol=search.run();
        
        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
       
    }
    if(!strcmp(argv[1],"mfo"))
    {
        // cout<<"Number of population:"<<argv[4]<<endl;
        mfo search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
        mfo::solution sol=search.run();
        
        // for(int j=0;j<sol.size();j++)
        //     cout<<sol[j]<<" ";
        // cout<<endl;
    }
    if(!strcmp(argv[1],"pso"))
    {
        cout<<"Number of population:"<<argv[4]<<endl;
        pso search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
        pso::solution sol=search.run();
        
        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
    }
    if(!strcmp(argv[1],"ga"))
    {
        cout<<"Number of population:"<<argv[4]<<endl;
        ga search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atof(argv[6]),atof(argv[7]),argv[8]);
        ga::dd1 sol=search.run();
        
        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
    }
    if(!strcmp(argv[1],"se"))
    {
        cout<<"Number of dimension:"<<argv[4]<<endl;
        cout<<"Number of regions:"<<argv[5]<<endl;
        cout<<"Number of searcher:"<<argv[6]<<endl;
        cout<<"Number of sample:"<<argv[7]<<endl;
        cout<<"Tournament Size:"<<argv[8]<<endl;
        
        se search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]));
        se::dd1 sol=search.run();
        
        for(int j=0;j<sol.size();j++)
            cout<<sol[j]<<" ";
        cout<<endl;
    }
    if(!strcmp(argv[1],"de"))
    {
        // cout<<"Number of population:"<<argv[4]<<endl;
        mfo search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
        mfo::solution sol=search.run();
        
        // for(int j=0;j<sol.size();j++)
        //     cout<<sol[j]<<" ";
        // cout<<endl;
    }
    if(!strcmp(argv[1],"mfode"))
    {
        // cout<<"Number of population:"<<argv[4]<<endl;
        mfo search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
        mfo::solution sol=search.run();
        
        // for(int j=0;j<sol.size();j++)
        //     cout<<sol[j]<<" ";
        // cout<<endl;
    }
    return 0;
}