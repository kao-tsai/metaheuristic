
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<string.h>
#include"hc.h"
using namespace std;

int main(int argc,char **argv){
    clock_t start,finish;
    double duration;
    cout<<"Algorithm:"<<argv[1]<<endl;
    cout<<"Number of runs:"<<argv[2]<<endl;
    cout<<"Number of iterations:"<<argv[3]<<endl;
    cout<<"Number of patterns:"<<argv[4]<<endl;
    cout<<"Name of file:"<<argv[5]<<endl;
    if(!strcmp(argv[1],"hc"))
    {
        
        hc search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5]);
        start=clock();
        hc::solution sol=search.run();
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        cout<<"final result:"<<endl;
        for(int i=0;i<sol.size();i++)
            cout<<sol[i];
        cout<<endl;
        cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    }

        return 0;
}
