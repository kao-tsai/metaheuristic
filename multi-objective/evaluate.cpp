#include<iostream>
#include<vector>
#include <string>
#include<fstream>
#include<cmath>
using namespace std;
double dis_cal(double x1,double y1,double x2,double y2);
void old_IGD(string function,int obj,int pop){
    ifstream opt;
    // FON_opt.txt
    opt.open("opt/"+function+"_opt.txt");
    
    vector<vector<double>> true_opt=vector<vector<double>> (1000,vector<double>(obj,0));
    double get_value;
    int i=0;
    while(i<pop){
        opt>>get_value;
        true_opt[i][0]=get_value;
        opt>>get_value;
        true_opt[i][1]=get_value;
        i++;
    }
    opt.close();
    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+function+"_mnsga2_0.txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    i=0;
    while(i<pop){
        opt>>get_value;
        all_sol[i][0]=get_value;
        opt>>get_value;
        all_sol[i][1]=get_value;
        i++;
    }
    open_sol.close();
    double min,temp;
    double sum=0;
    for(i=0;i<pop;i++)
    {
        min=dis_cal(all_sol[i][0],all_sol[i][1],true_opt[0][0],true_opt[0][1]);
        for(int j=1;j<1000;j++)
        {
            temp=dis_cal(all_sol[i][0],all_sol[i][1],true_opt[j][0],true_opt[j][1]);
            if(temp<min){
                min=temp;
            }
        }
        sum+=min;
    }
    
    cout<<"Function("<<function<<") IGD="<<sum/pop<<endl;
}

double dis_cal(double x1,double y1,double x2,double y2){
    return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}