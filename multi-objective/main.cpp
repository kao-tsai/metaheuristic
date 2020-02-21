#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include"evaluate.cpp"
// #include"MOPSO/mopso.h"
#include"MOPSO/mopso(noWH).h"
using namespace std;
/*
vector<int>test1(vector<int> th){
    th[0]=50;
    th[1]=49;
    th[2]=48;
    th.push_back(47);
    return th;
}*/
int main(int argc,char **argv)
{
    clock_t start,finish;
    double duration;
    cout<<"Algorithm:"<<argv[1]<<endl;                      //nsga2
    cout<<"Number of runs:"<<argv[2]<<endl;          //default 30
    cout<<"Number of iterations:"<<argv[3]<<endl;//default 250
    cout<<"Number of population:"<<argv[4]<<endl;                  //default 500
    
    if(!strcmp(argv[1],"bnsga2"))
    {
        // bnsgaii search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atof(argv[6]),argv[7]);
        // vector<vector<double>> sol=search.run();
        /*
        for(int i=0;i<sol.size();i++)
            cout<<sol[i];
        cout<<endl;*/
        // for(int i=0;i<sol.size();i++){
        //     for(int j=0;j<sol[i].size();j++)
        //         cout<<sol[i][j]<<" ";
        //     cout<<endl;
        // }
        old_IGD("FON",2,500);
        old_IGD("SCH",2,500);
        old_IGD("ZDT1",2,500);
        old_IGD("ZDT2",2,500);
        old_IGD("ZDT3",2,500);
        old_IGD("ZDT4",2,100);
        old_IGD("ZDT6",2,500);
    }
    else if(!strcmp(argv[1],"mopso")){

        start=clock();
        mopso search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atoi(argv[6]),argv[7]);
        mopso::population sol=search.run();
        finish=clock();
        duration=(double)(finish-start)/CLOCKS_PER_SEC;
        
        for(int i=0;i<sol.size();i++){
            for(int j=0;j<sol[i].size();j++)
                cout<<sol[i][j]<<" ";
            cout<<endl;
        }
        cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    }

    return 0;
}
