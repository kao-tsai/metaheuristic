#include<iostream>
#include<vector>
#include <string>
#include<fstream>
#include<cmath>
#include<limits>
using namespace std;

void quick(vector<double>&,vector<int>&,int,int);

double dis_cal(vector<double> sol1,vector<double>sol2);
double d_plus_cal(vector<double>sol,vector<double>true_sol);
double NSGAII_GD(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
    ifstream opt;
    //opt/和nsgaiiTest/要互換
    opt.open(opt_data+"/"+function+"_opt.txt");
    // opt.open("nsgaiiTest/"+function+"_opt.txt");
    
    vector<vector<double>> true_opt=vector<vector<double>> (opt_num,vector<double>(obj,0));
    double get_value;
    int i=0;
    while(i<opt_num){
        opt>>get_value;
        true_opt[i][0]=get_value;
        // cout<<get_value<<" ";
        opt>>get_value;
        true_opt[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    opt.close();
    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    i=0;
    while(i<pop){
        open_sol>>get_value;
        all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        open_sol>>get_value;
        all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    open_sol.close();
    double min,temp;
    double sum=0;
    for(i=0;i<pop;i++)
    {
        min=dis_cal(all_sol[i],true_opt[0]);
        for(int j=1;j<opt_num;j++)
        {
            temp=dis_cal(all_sol[i],true_opt[j]);
            if(temp<min){
                min=temp;
            }
        }
        sum+=min;
    }
    // ofstream output_IGD;
    // output_IGD.open("pareto/"+function+"/nsgaii/"+function+"_IGD");
    // output_IGD<<sum/pop<<endl;
    // output_IGD.close();
    // cout<<"Function("<<function<<") IGD="<<sum/pop<<endl;
    return sum/(double)pop;
}
//Generational Distance
double GD(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
    ifstream opt;
    //opt/和nsgaiiTest/要互換
    opt.open(opt_data+"/"+function+"_opt.txt");
    // opt.open("nsgaiiTest/"+function+"_opt.txt");
    
    vector<vector<double>> true_opt=vector<vector<double>> (opt_num,vector<double>(obj,0));
    double get_value;
    int i=0;
    while(i<opt_num){
        opt>>get_value;
        true_opt[i][0]=get_value;
        // cout<<get_value<<" ";
        opt>>get_value;
        true_opt[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    opt.close();
    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    i=0;
    while(i<pop){
        open_sol>>get_value;
        all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        open_sol>>get_value;
        all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    open_sol.close();
    double min,temp;
    double sum=0;
    for(i=0;i<pop;i++)
    {
        min=dis_cal(all_sol[i],true_opt[0]);
        for(int j=1;j<opt_num;j++)
        {
            temp=dis_cal(all_sol[i],true_opt[j]);
            if(temp<min){
                min=temp;
            }
        }
        sum+=pow(min,2.0);
    }
    // ofstream output_IGD;
    // output_IGD.open("pareto/"+function+"/nsgaii/"+function+"_IGD");
    // output_IGD<<sum/pop<<endl;
    // output_IGD.close();
    // cout<<"Function("<<function<<") IGD="<<sum/pop<<endl;
    return sqrt(sum)/(double)pop;
}
//Inverted Generational Distance
double IGD(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
    ifstream opt;
    //opt/和nsgaiiTest/要互換
    opt.open(opt_data+"/"+function+"_opt.txt");
    // opt.open("nsgaiiTest/"+function+"_opt.txt");
    
    vector<vector<double>> true_opt=vector<vector<double>> (opt_num,vector<double>(obj,0));
    double get_value;
    int i=0;
    // while(i<opt_num){
    //     opt>>get_value;
    //     true_opt[i][0]=get_value;
    //     // cout<<get_value<<" ";
    //     opt>>get_value;
    //     true_opt[i][1]=get_value;
    //     // cout<<get_value<<endl;
    //     i++;
    // }
    opt_num=0;
    while(opt>>get_value){
        true_opt[opt_num][0]=get_value;
        for(int k=1;k<obj;k++){
            opt>>get_value;
            true_opt[opt_num][k]=get_value;
        }
        opt_num++;

    }
    opt.close();

    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    i=0;
    while(i<pop){
        for(int j=0;j<obj;j++){
            open_sol>>get_value;
            all_sol[i][j]=get_value;
        }
        // open_sol>>get_value;
        // all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        // open_sol>>get_value;
        // all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    open_sol.close();
    double min,temp;
    double p=1;
    double sum=0;
    for(i=0;i<opt_num;i++)
    {
        min=dis_cal(all_sol[0],true_opt[i]);
        for(int j=1;j<pop;j++)
        {
            temp=dis_cal(all_sol[j],true_opt[i]);
            if(temp<min)
                min=temp;
                                                                                                             
        }
        sum+=pow(min,p);
    }
    // ofstream output_IGD;
    // output_IGD.open("pareto/"+function+"/nsgaii/"+function+"_IGD");
    // output_IGD<<sum/pop<<endl;
    // output_IGD.close();
    // cout<<"Function("<<function<<") IGD="<<sum/pop<<endl;
    return pow(sum,(1.0/p))/(double)opt_num;
}
double IGD_plus(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
    ifstream opt;
    //opt/和nsgaiiTest/要互換
    opt.open(opt_data+"/"+function+"_opt.txt");
    
    vector<vector<double>> true_opt=vector<vector<double>> (opt_num,vector<double>(obj,0));
    double get_value;
    int i=0;
    opt_num=0;
    while(opt>>get_value){
        true_opt[opt_num][0]=get_value;
        for(int k=1;k<obj;k++){
            opt>>get_value;
            true_opt[opt_num][k]=get_value;
        }
        opt_num++;

    }
    opt.close();
    // cout<<opt_num<<endl;
    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    i=0;
    while(i<pop){
        open_sol>>get_value;
        all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        open_sol>>get_value;
        all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    open_sol.close();
    double min,temp;
    double p=2;
    double sum=0;
    for(i=0;i<opt_num;i++)
    {
        min=d_plus_cal(all_sol[0],true_opt[i]);
        for(int j=1;j<pop;j++)
        {
            temp=d_plus_cal(all_sol[j],true_opt[i]);
            if(temp<min)
                min=temp;
        }
        sum+=min;
    }
    // ofstream output_IGD;
    // output_IGD.open("pareto/"+function+"/nsgaii/"+function+"_IGD");
    // output_IGD<<sum/pop<<endl;
    // output_IGD.close();
    // cout<<"Function("<<function<<") IGD="<<sum/pop<<endl;
    return pow((sum/(double)opt_num),1.0/p);
}
//以下SP衡量方式取自NSGAII，若換演算法則須修改
double NSGAII_SP(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
    ifstream opt;
    // FON_opt.txt
    opt.open(opt_data+"/"+function+"_opt.txt");
    // opt.open("nsgaiiTest/"+function+"_opt.txt");
    vector<double> min_obj(obj,numeric_limits<double>::infinity());
    vector<vector<double>> extreme_sol(obj,vector<double>());
    vector<vector<double>> true_opt=vector<vector<double>> (opt_num,vector<double>(obj,0));
    double get_value;
    int i=0;
    while(i<opt_num){
        opt>>get_value;
        true_opt[i][0]=get_value;
 
        opt>>get_value;
        true_opt[i][1]=get_value;

        if(min_obj[0]>true_opt[i][0]){
            min_obj[0]=true_opt[i][0];
            extreme_sol[0]=true_opt[i];
        }
        if(min_obj[1]>true_opt[i][1]){
            min_obj[1]=true_opt[i][1];
            extreme_sol[1]=true_opt[i];
        }
        i++;
    }
    opt.close();

    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    vector<double> obj_one(pop);
    vector<int> sort_pos(pop);
    i=0;
    while(i<pop){
        open_sol>>get_value;
        all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        open_sol>>get_value;
        all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        obj_one[i]=all_sol[i][0];
        //
        // cout<<obj_one[i]<<endl;

        sort_pos[i]=i;
        i++;
    }
    open_sol.close();

    quick(obj_one,sort_pos,0,pop-1);
    // cout<<obj_one[0]<<endl;
    // cout<<sort_pos[0]<<endl;
    // exit(1);

    vector<double> d_n(pop+1);
    //計算所有相鄰解的距離以及邊界解的距離
    d_n[0]=dis_cal(extreme_sol[0],all_sol[sort_pos[0]]);
    i=1;
    for(int j=0;j<pop-1;j++){
        d_n[i]=dis_cal(all_sol[sort_pos[j]],all_sol[sort_pos[j+1]]);
        i++;
    }
    d_n[pop]=dis_cal(extreme_sol[1],all_sol[sort_pos[pop-1]]);
    //計算出平均距離
    double average_d=0;
    for(int j=1;j<pop;j++)
        average_d+=d_n[j];
    average_d/=(double)(pop-1);
    //開始計算Spacing
    double sp_value=0;
    for(int j=1;j<pop;j++)
        sp_value+=fabs(d_n[j]-average_d);
    sp_value=(sp_value+d_n[0]+d_n[pop])/(d_n[0]+d_n[pop]+(pop-1)*average_d);

    return sp_value;
}
//以下SP衡量方式取自
double SP(string function,string algo,string opt_data,int run,int obj,int pop,int opt_num){
   
    double get_value;

    ifstream open_sol;
    open_sol.open("pareto/"+function+"/"+algo+"/"+function+"_"+algo+"_"+to_string(run)+".txt");
    vector<vector<double>> all_sol=vector<vector<double>> (pop,vector<double>(obj,0));
    // vector<double> obj_one(pop);
    vector<int> sort_pos(pop);
    int i=0;
    while(i<pop){
        for(int j=0;j<obj;j++){
            open_sol>>get_value;
            all_sol[i][j]=get_value;
        }
        // open_sol>>get_value;
        // all_sol[i][0]=get_value;
        // cout<<get_value<<" ";
        // open_sol>>get_value;
        // all_sol[i][1]=get_value;
        // cout<<get_value<<endl;
        i++;
    }
    open_sol.close();

    // quick(obj_one,sort_pos,0,pop-1);
    

    // vector<double> d_n(pop+1);
    double min=std::numeric_limits<double>::infinity();
    vector<double> di(pop,min);
    double temp_sum=0;

    double average_sum=0;
    //計算所有相鄰解的距離以及邊界解的距離
    for(i=0;i<pop;i++){
        for(int j=0;j<pop;j++){
            if(i!=j){
                temp_sum=0;
                for(int k=0;k<obj;k++)
                    temp_sum+=fabs(all_sol[i][k]-all_sol[j][k]);
                if(temp_sum<di[i])
                    di[i]=temp_sum;
            }
        }
        average_sum+=di[i];
    }
    //計算出平均距離
    double average_d=average_sum/(double)pop;
   
    //開始計算Spacing
    double sp_value=0;
    for(int i=0;i<pop;i++)
        sp_value+=pow(di[i]-average_d,2.0);
    sp_value/=pop;

    return sqrt(sp_value);
}

void quick(vector<double>& obj_one,vector<int>& sort_pos,int left,int right){
    int index;
    double double_tmp;
    int int_tmp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rand()%(right-left+1)+left;

        double_tmp = obj_one[right];
        obj_one[right] = obj_one[index];
        obj_one[index] = double_tmp;

        int_tmp = sort_pos[right];
        sort_pos[right] = sort_pos[index];
        sort_pos[index] = int_tmp;

        pivot = obj_one[right];

        i = left-1;
        for (j=left; j<right; j++)
        {
            if (obj_one[j] <= pivot)
            {
                i+=1;
                double_tmp = obj_one[j];
                obj_one[j] = obj_one[i];
                obj_one[i] = double_tmp;
                
                int_tmp = sort_pos[j];
                sort_pos[j] = sort_pos[i];
                sort_pos[i] = int_tmp;
            }
        }
        index=i+1;
        double_tmp = obj_one[index];
        obj_one[index] = obj_one[right];
        obj_one[right] = double_tmp;

        int_tmp = sort_pos[index];
        sort_pos[index] = sort_pos[right];
        sort_pos[right] = int_tmp;

        quick(obj_one,sort_pos, left, index-1);
        quick(obj_one,sort_pos, index+1, right);
    }
}

double dis_cal(vector<double> sol1,vector<double>sol2){
    double sum=0.0;
    int sol_size=sol1.size();
    for(int i=0;i<sol_size;i++){
        sum+=pow(sol1[i]-sol2[i],2.0);
    }
    return sqrt(sum);
}
double d_plus_cal(vector<double>sol,vector<double>true_sol){
    double sum=0,temp;
    for(int i=0;i<sol.size();i++){
        temp=sol[i]-true_sol[i];
        if(temp>0)
            sum+=pow(temp,2.0);
        else 
            sum+=0;
    }
    return sqrt(sum);
}