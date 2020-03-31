#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include"evaluate.cpp"
// #include"NSGAII/bnsgaii.h"
#include"SE/se.h"
// #include"MOPSO/mopso.h"
// #include"MOPSO/mopso(noWH).h"
#include<fstream>
#include <math.h>
#include <algorithm>
using namespace std;

typedef vector<int> d1;
typedef vector<d1> d2;
typedef vector<double> dd1;
typedef vector<dd1> dd2;
dd2 combine_fit(103,dd1(2));
void quick_sort(dd1& value,d1& this_rank,int left,int right){
    int index;
    int int_tmp;
    double double_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        double_tmp = value[right];
        value[right] = value[index];
        value[index] = double_tmp;

        int_tmp = this_rank[right];
        this_rank[right] = this_rank[index];
        this_rank[index] = int_tmp;

        pivot = value[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (value[j] <= pivot)
            {
                i+=1;
                double_tmp = value[j];
                value[j] = value[i];
                value[i] = double_tmp;
                
                int_tmp = this_rank[j];
                this_rank[j] = this_rank[i];
                this_rank[i] = int_tmp;
            }
        }
        index=i+1;
        double_tmp = value[index];
        value[index] = value[right];
        value[right] = double_tmp;

        int_tmp = this_rank[index];
        this_rank[index] = this_rank[right];
        this_rank[right] = int_tmp;

        quick_sort(value,this_rank, left, index-1);
        quick_sort(value,this_rank, index+1, right);
    }
}
void quick_sort_a(d1& value,dd1& this_rank,int left,int right){
    int index;
    int int_tmp;
    double double_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        int_tmp = value[right];
        value[right] = value[index];
        value[index] = double_tmp;

        double_tmp = this_rank[right];
        this_rank[right] = this_rank[index];
        this_rank[index] = int_tmp;

        pivot = value[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (value[j] <= pivot)
            {
                i+=1;

                int_tmp = value[j];
                value[j] = value[i];
                value[i] = double_tmp;

                double_tmp = this_rank[j];
                this_rank[j] = this_rank[i];
                this_rank[i] = int_tmp;
                
            }
        }
        index=i+1;
        int_tmp = value[index];
        value[index] = value[right];
        value[right] = double_tmp;

        double_tmp = this_rank[index];
        this_rank[index] = this_rank[right];
        this_rank[right] = int_tmp;

        quick_sort_a(value,this_rank, left, index-1);
        quick_sort_a(value,this_rank, index+1, right);
    }
}
void quick_sort_obj(int obj_level,d1& obj_id,d1& sorting_id,int left,int right){
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        temp = obj_id[right];
        obj_id[right] = obj_id[index];
        obj_id[index] = temp;

        temp = sorting_id[right];
        sorting_id[right] = sorting_id[index];
        sorting_id[index] = temp;

        pivot = combine_fit[obj_id[right]][obj_level];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (combine_fit[obj_id[j]][obj_level] <= pivot)
            {
                i+=1;
                temp = obj_id[j];
                obj_id[j] = obj_id[i];
                obj_id[i] = temp;
                
                temp = sorting_id[j];
                sorting_id[j] = sorting_id[i];
                sorting_id[i] = temp;
            }
        }
        index=i+1;
        temp = obj_id[index];
        obj_id[index] = obj_id[right];
        obj_id[right] = temp;

        temp = sorting_id[index];
        sorting_id[index] = sorting_id[right];
        sorting_id[right] = temp;

        quick_sort_obj(obj_level,obj_id,sorting_id, left, index-1);
        quick_sort_obj(obj_level,obj_id,sorting_id, index+1, right);
    }
}
void cda(d1& this_rank){
    int front_size=this_rank.size();
    vector<double> init_dis(front_size,0.0);

    d2 sorting_id(2,d1(front_size));
    d2 obj_id(2);
    
    for(int i=0;i<2;i++){
        for(int j=0;j<front_size;j++)
            sorting_id[i][j]=j;
        obj_id[i]=this_rank;
        quick_sort_obj(i,obj_id[i],sorting_id[i],0,front_size-1);
    }

    //計算擁擠程度
    for (int i=0; i<2; i++)
        init_dis[sorting_id[i][0]]= std::numeric_limits<double>::infinity();
    
    for(int i=0;i<2;i++)
        for(int j=1;j<front_size-1;j++){
            if(init_dis[j]!=numeric_limits<double>::infinity()){
                if(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]== combine_fit[this_rank[sorting_id[i][0]]][i])
                    init_dis[sorting_id[i][j]]+=0.0;
                else
                    init_dis[sorting_id[i][j]]+=(combine_fit[this_rank[sorting_id[i][j+1]]][i]-combine_fit[this_rank[sorting_id[i][j-1]]][i])/(combine_fit[this_rank[sorting_id[i][front_size-1]]][i]-combine_fit[this_rank[sorting_id[i][0]]][i]);
            }
        }
    //距離密度值進行正規化
    for (int j=0; j<front_size; j++){
        if (init_dis[j]!= numeric_limits<double>::infinity())
           init_dis[j] = init_dis[j]/2.0;
        //    cout<<"init_dis "<<j<<": "<<init_dis[j]<<endl;
    }
    //
    d1 tmp_rank=this_rank;

    // quick_sort(init_dis,this_rank,0,front_size-1);
    // reverse(init_dis.begin(),init_dis.end());
    // reverse(this_rank.begin(),this_rank.end());
    // quick_sort_a(this_rank,init_dis,0,front_size-1);
    for(int i=0;i<front_size;i++){
        cout<<this_rank[i]<<",cd= "<<init_dis[i]<<endl;
    }
}

int main(int argc,char **argv)
{
    clock_t start,finish;
    double duration;
    // cout<<"Algorithm:"<<argv[1]<<endl;                      //nsga2
    // cout<<"Number of runs:"<<argv[2]<<endl;          //default 30
    // cout<<"Number of iterations:"<<argv[3]<<endl;//default 250
                 
   
    if(!strcmp(argv[1],"bnsga2"))
    {
        // cout<<"Number of population:"<<argv[4]<<endl;//default 500    
        
        // bnsgaii search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atof(argv[6]),argv[7]);
        
        // search.run();

        // int run_num=atoi(argv[2]);
        // int pop_num=atoi(argv[4]);
        // string func_name=argv[7];
        // double tmp;
        // double sum_IGD=0;
        // double sum_SP=0;
        // double sum_MS=0;
        // ofstream output_igd;
        // output_igd.open("pareto/"+func_name+"/bnsga2/IGD.txt");
        // for(int i=0;i<run_num;i++){
        //     tmp=IGD(func_name,"bnsga2",i,2,pop_num,500);
        //     sum_IGD+=tmp;
        //     output_igd<<tmp<<endl;
        // }
        // output_igd<<sum_IGD/run_num<<endl;
        // output_igd.close();

        // ofstream output_sp;
        // output_sp.open("pareto/"+func_name+"/bnsga2/SP.txt");
        // for(int i=0;i<run_num;i++){
        //     tmp=SP(func_name,"bnsga2",i,2,pop_num,500);
        //     sum_SP+=tmp;
        //     output_sp<<tmp<<endl;
        // }
        // output_sp<<sum_SP/run_num<<endl;
        // output_sp.close();
        
    }

    // else if(!strcmp(argv[1],"mopso")){
    //     cout<<"Number of population:"<<argv[4]<<endl;     
    //     start=clock();
    //     mopso search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atof(argv[5]),atoi(argv[6]),argv[7]);
    //     mopso::population sol=search.run();
    //     finish=clock();
    //     duration=(double)(finish-start)/CLOCKS_PER_SEC;
        
    //     for(int i=0;i<sol.size();i++){
    //         for(int j=0;j<sol[i].size();j++)
    //             cout<<sol[i][j]<<" ";
    //         cout<<endl;
    //     }
    //     cout<<fixed<<setprecision(4)<<duration<<".sec"<<endl;
    // }
     if(!strcmp(argv[1],"se"))
    {
        
        se search(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),argv[9],atof(argv[10]),atof(argv[11]));
        search.run();
        
        // for(int i=0;i<sol.size();i++){
        //     for(int j=0;j<sol[i].size();j++)
        //         cout<<sol[i][j]<<" ";
        //     cout<<endl;
        // }
        
        int run_num=atoi(argv[2]);
        int pop_num=atoi(argv[4]);
        string func_name=argv[9];
        double tmp;
        double sum_IGD=0;
        double sum_SP=0;
        double sum_MS=0;
        ofstream output_igd;
        output_igd.open("pareto/"+func_name+"/se/IGD.txt");
        for(int i=0;i<run_num;i++){
            tmp=IGD(func_name,"se",i,2,pop_num,500);
            sum_IGD+=tmp;
            output_igd<<tmp<<endl;
        }
        output_igd<<sum_IGD/run_num<<endl;
        output_igd.close();

        ofstream output_sp;
        output_sp.open("pareto/"+func_name+"/se/SP.txt");
        for(int i=0;i<run_num;i++){
            tmp=SP(func_name,"se",i,2,pop_num,500);
            sum_SP+=tmp;
            output_sp<<tmp<<endl;
        }
        output_sp<<sum_SP/run_num<<endl;
        output_sp.close();
        
        // vector<vector<double>> fit={{0.0395688,0.961342},{0.0400102,0.96118}};
        // double middle=0.0;
        // //計算擁擠程度
        
        // middle+=(fit[1][0]-fit[0][0])/(0.981717-1.20078e-06);
        // middle+=(fit[0][1]-fit[1][1])/(0.981749-6.00214e-07);
        // cout<<middle/2<<endl;
        // ifstream test_6_in;
        // test_6_in.open("test_delete/FON_se_full_29.txt");
        // d1 rank(103);
        // for(int i=0;i<103;i++){
        //     test_6_in>>combine_fit[i][0];
        //     test_6_in>>combine_fit[i][1];
        //     rank[i]=i;
        // }
        // test_6_in.close();
        // cda(rank);

    }
    
    return 0;

}
