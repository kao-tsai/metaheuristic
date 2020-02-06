#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>

#include"NSGAII/bnsgaii.h"
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
    cout<<"Algorithm:"<<argv[1]<<endl;                      //nsga2
    cout<<"Number of runs:"<<argv[2]<<endl;          //default 30
    cout<<"Number of iterations:"<<argv[3]<<endl;//default 250
    cout<<"Number of population:"<<argv[4]<<endl;                  //default 500
    
    if(!strcmp(argv[1],"bnsga2"))
    {
        bnsgaii search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atof(argv[6]),argv[7]);
        vector<vector<double>> sol=search.run();
        /*
        for(int i=0;i<sol.size();i++)
            cout<<sol[i];
        cout<<endl;*/
        for(int i=0;i<sol.size();i++){
            for(int j=0;j<sol[i].size();j++)
                cout<<sol[i][j]<<" ";
            cout<<endl;
        }

        /*
        typedef vector<int> solution;
        typedef vector<solution> population;
        solution test(3);
        population two;
        test[0]=1;
        test[1]=1;
        test[2]=1;
        test=test1(test);
        cout<<test[0]<<endl;
        cout<<test[1]<<endl;
        cout<<test[2]<<endl;
        cout<<test[3]<<endl;
        two.push_back(test);
        two.push_back(test);
        cout<<two.size()<<endl;
        two[0].insert(two[0].begin(),test.begin(),test.end());
        cout<<two[0].size()<<endl;
        population two2(2);
        two.insert(two.end(),two2.begin(),two2.end());
        cout<<two.size()<<endl;
        two=two2;
        cout<<two.size()<<endl;*/
        
    }

    return 0;
}
