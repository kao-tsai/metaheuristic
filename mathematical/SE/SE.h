//Based on HHO's transition
//Divide the regions based on Dimension----------> Bound:0~100:R1=[(0,50),(0,50)],R2=[(50,100),(0,50)]
//                                                                                                                                R3=[(0,50),(50,100)],R4=[(50,100),(50,100)]

#ifndef __SE_H_INCLUDED__
#define __SE_H_INCLUDED__
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <math.h>
#include<iomanip>
#include<random>
#include<algorithm>
#include<functional>
using namespace std;

class se{
	public:
		typedef vector<int> d1;
        typedef vector<d1> d2;
        typedef vector<d2> d3;
        typedef vector<double> dd1;
        typedef vector<dd1> dd2;
        typedef vector<dd2> dd3;

		se(int,int,int,int,int,int,int);
		dd1 run();
	private:
		void init();
        double bound_rand(double,double);
		double evaluation(dd1&);
        void divide_region();
        void arrange();
        void cal_ev();
        void choose_sample();
        void select_player();
        void HHO_transit(double);
        void marketing_survey();
        dd1 levy_flight();
        dd1 cal_mean(dd2&);
        dd1 cal_mean(dd3&);

private:
	int numRuns;
	int numIter;
	int numPatterns;
    
    dd2 searcher;
    dd3 sample;
    vector<dd3> sampleV;
    dd3 region_bound;
    dd1 ta;
    dd1 tb;
    d1 searcher_belong;
    dd2 EV;
    dd2 region_best;
    dd1 region_best_fit;
    dd2 sample_fit;
    dd3 sampleV_fit;
    dd1 searcher_fit;
    dd1 best_sol;
    dd1 upperbound;
    dd1 lowerbound;

    dd1 gnu_eva;

    double inf;
    int clip_bit_num;
    double best_obj_val;
    int region_num;
    int searcher_num;
    int sample_num;
    int tour_size;
    int evaluate_num;
    int not_evaluate;
    int max_evaluate;
};

se::se(int xNumRuns,
int xNumIter,
int xdim,
int xRegion,
int in_searcher_num,
int in_sample_num,
int in_tournament_size
)
{
    srand(time(0));
    inf=std::numeric_limits<double>::infinity();
    numRuns=xNumRuns;
	numIter=xNumIter;
	numPatterns=xdim;
    region_num=xRegion;
    searcher_num=in_searcher_num;
    sample_num=in_sample_num;
    tour_size=in_tournament_size;
    max_evaluate=270000;
}


void se::init(){
    best_obj_val=inf;
    
    evaluate_num=0;
    gnu_eva[evaluate_num]+=best_obj_val;
    best_sol.assign(numPatterns,0);
    searcher_belong.assign(searcher_num,0);
    searcher.assign(searcher_num,dd1(numPatterns,0));
    sample.assign(region_num,dd2(sample_num,dd1(numPatterns,0)));
    sampleV.assign(searcher_num, dd3(region_num, dd2(sample_num, dd1(numPatterns, 0))));
    sample_fit.assign(region_num, dd1(sample_num, 0.0));
    sampleV_fit.assign(searcher_num,dd2(region_num,dd1(sample_num,0.0)));
    searcher_fit.assign(searcher_num,0.0);
    region_best_fit.assign(region_num,inf);
    region_best.assign(region_num,dd1(numPatterns,0.0));
    //自行設定特定公式的邊界
    upperbound.assign(numPatterns,30);
    lowerbound.assign(numPatterns,-30);
	for (int i = 0; i < searcher_num; i++)
        for(int j=0;j<numPatterns;j++)
            searcher[i][j]=bound_rand(lowerbound[j],upperbound[j]);
}
void se::divide_region(){
    //Lower Bound在上，Upper Bound在下
    region_bound.assign(region_num,dd2(2,dd1(clip_bit_num,0)));
    for(int i=0;i<region_num;i++){
        int bin=i;     
        for(int j=clip_bit_num-1;j>=0;j--){
            double divide_two=(upperbound[j]+lowerbound[j])/2;
            if(bin%2==0){
                region_bound[i][0][j]=lowerbound[j];
                region_bound[i][1][j]=divide_two+(upperbound[j]-lowerbound[j]);//10
            }
            else if(bin%2==1){
                region_bound[i][0][j]=divide_two-(upperbound[j]-lowerbound[j]);//10
                region_bound[i][1][j]=upperbound[j];
            }
            bin/=2;
        }
    }
    //記得刪-------------------------------------------------------------------------------------//
    // for(int i=0;i<region_num;i++)
    // {
    //     for(int j=0;j<clip_bit_num;j++)
    //     {
    //         cout<<"Region "<<i<<" lowerbound:"<<region_bound[i][0][j]<<endl;
    //         cout<<"Region "<<i<<" upperbound:"<<region_bound[i][1][j]<<endl;
    //         cout<<endl;
    //     }
    // }
    // exit(1);
    //------------------------------------------------------------------------------------------//
}
double se::bound_rand(double lw,double ub){
    return ((double)rand()/RAND_MAX)*(ub-lw)+lw;
}
void se::arrange(){
    clip_bit_num=log2(region_num);

    divide_region();
    for(int i=0;i<searcher_num;i++){
        searcher_belong[i]=i%region_num;
        for(int j=0;j<clip_bit_num;j++)
            searcher[i][j]=bound_rand(region_bound[searcher_belong[i]][0][j],region_bound[searcher_belong[i]][1][j]);
        for(int j=clip_bit_num;j<numPatterns;j++)
            searcher[i][j]=bound_rand(lowerbound[j],upperbound[j]);
    }

    //Sample初始化------------------->隨機Sample
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++)
            for(int k=0;k<numPatterns;k++)
                if(k<clip_bit_num)
                    sample[i][j][k]=bound_rand(region_bound[i][0][k],region_bound[i][1][k]);
                else
                    sample[i][j][k]=bound_rand(lowerbound[k],upperbound[k]);
    
    //SampleV初始化----------------->將所有Sample複製進來
    for(int i=0;i<searcher_num;i++)
        sampleV[i]=sample;
    
    tb.assign(region_num,0.0);
    ta.assign(region_num,1.0);
    for(int i=0;i<searcher_num;i++){
        tb[searcher_belong[i]]++;
        ta[searcher_belong[i]]=1.0;
    }
                
    EV.assign(searcher_num,dd1(region_num,0.0));
    //
    //
    //
}
//幾種問題：1.考量是否由Rabbit當Searcher  2.質心由(1)該區域計算為主或是(2)全區域總和的平均
//選擇方法：質心暫時只選該區域

void se::HHO_transit(double E1){
    double E0,escaping_energy,q;
    
    for(int i=0;i<searcher_num;i++){
        for(int j=0;j<region_num;j++){
            for(int k=0;k<sample_num;k++){
                E0=2.0*((double)rand()/RAND_MAX)-1;
                escaping_energy=E0*E1;
                //記住刪除-------------//
                // cout<<"escaping_energy="<<escaping_energy<<endl;
                //------------------------//
                if(fabs(escaping_energy)>=1){
                    q=(double)rand()/RAND_MAX;
                    if(q<0.5){
                        dd1 rand_hawk=sampleV[i][j][rand()%sample_num];
                        double r1=(double)rand()/RAND_MAX;
                        double r2=(double)rand()/RAND_MAX;
                        for(int l=0;l<numPatterns;l++)
                            sampleV[i][j][k][l]=rand_hawk[l]-r1*fabs(rand_hawk[l]-2.0*r2*sampleV[i][j][k][l]); 
                    }
                    else if(q>=0.5){
                        double r1=(double)rand()/RAND_MAX;
                        double r2=(double)rand()/RAND_MAX;
                        dd1 mean=cal_mean(sampleV[i][j]);
                        // dd1 mean=cal_mean(sampleV[i]);
                        for(int l=0;l<numPatterns;l++)
                            if(l<clip_bit_num)
                                sampleV[i][j][k][l]=(searcher[i][l]-mean[l])-r1*((region_bound[j][1][l]-region_bound[j][0][l])*r2+region_bound[j][0][l]);
                            else
                                sampleV[i][j][k][l]=(searcher[i][l]-mean[l])-r1*((upperbound[l]-lowerbound[l])*r2+lowerbound[l]);
                    }
                }
                else if(fabs(escaping_energy)<1){
                    //do the exploitation phase
                    
                    double r=(double)rand()/RAND_MAX;
                    if(r>=0.5 && fabs(escaping_energy)<0.5){
                        for(int l=0;l<numPatterns;l++)
                            sampleV[i][j][k][l]=searcher[i][l]-escaping_energy*fabs(searcher[i][l]-sampleV[i][j][k][l]);
                    }
                    if(r>=0.5 && fabs(escaping_energy)>=0.5){
                        double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                        for(int l=0;l<numPatterns;l++)
                            sampleV[i][j][k][l]=(searcher[i][l]-sampleV[i][j][k][l])-escaping_energy*fabs(jump_strength*searcher[i][l]-sampleV[i][j][k][l]);
                    }
                    if(r<0.5 && fabs(escaping_energy)>=0.5){
                        //
                        double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                        dd1 get_levy=levy_flight();
                        
                        dd1 X1(numPatterns),X2(numPatterns);
                        for(int l=0;l<numPatterns;l++)
                            X1[l]=searcher[i][l]-escaping_energy*fabs(jump_strength*searcher[i][l]-sampleV[i][j][k][l]);
                        if(evaluation(X1)<evaluation(sampleV[i][j][k]))
                            sampleV[i][j][k]=X1;
                        else
                        {
                            for(int l=0;l<numPatterns;l++)
                                X2[l]=X1[l]+((double)rand()/RAND_MAX)*get_levy[l];
                            if(evaluation(X2)<evaluation(sampleV[i][j][k]))
                                sampleV[i][j][k]=X2;
                        }

                    }
                    if(r<0.5 && fabs(escaping_energy)<0.5){
                        //
                        double jump_strength=2.0*(1-(double)rand()/RAND_MAX);
                        dd1 get_levy=levy_flight();
                        dd1 X1(numPatterns),X2(numPatterns);
                        dd1 mean=cal_mean(sampleV[i][j]);
                        // dd1 mean=cal_mean(sampleV[i]);
                        for(int l=0;l<numPatterns;l++)
                            X1[l]=searcher[i][l]-escaping_energy*fabs(jump_strength*searcher[i][l]-mean[l]);
                        if(evaluation(X1)<evaluation(sampleV[i][j][k]))
                            sampleV[i][j][k]=X1;
                        else
                        {
                            for(int l=0;l<numPatterns;l++)
                                X2[l]=X1[l]+((double)rand()/RAND_MAX)*get_levy[l];
                            if(evaluation(X2)<evaluation(sampleV[i][j][k]))
                                sampleV[i][j][k]=X2;
                        }
                    }
                }
                
                //檢查是否超出邊界或超出所屬範圍(Region)
                for(int l=0;l<numPatterns;l++)
                    if(l<clip_bit_num){
                        if(sampleV[i][j][k][l]<region_bound[j][0][l])
                            sampleV[i][j][k][l]=region_bound[j][0][l];
                        else if(sampleV[i][j][k][l]>region_bound[j][1][l])
                            sampleV[i][j][k][l]=region_bound[j][1][l];
                    }
                    else{
                        if(sampleV[i][j][k][l]<lowerbound[l])
                            sampleV[i][j][k][l]=lowerbound[l];
                        else if(sampleV[i][j][k][l]>upperbound[l])
                            sampleV[i][j][k][l]=upperbound[l];
                    }         
                
            }
        }
    }
}

se::dd1 se::levy_flight(){
    default_random_engine gen;
    normal_distribution<double> randn(0.0,1.0);
    double beta=1.5;
    double sigma=(gamma(1+beta)*sin((beta/2.0)*M_PI)/(gamma((1+beta)/2)*beta*pow((beta-1)/2,2)));
    sigma=pow(sigma,(1/beta));
    dd1 u(numPatterns);
    for(int i=0;i<numPatterns;i++)
        u[i]=(randn(gen)*sigma)/pow(fabs(randn(gen)),1.0/beta);
    return u;
}

se::dd1 se::cal_mean(dd2& all_sol){
    dd1 sum(numPatterns,0.0);
    int population_num=all_sol.size();
    for(int i=0;i<population_num;i++)
        transform(sum.begin(),sum.end(),all_sol[i].begin(),sum.begin(),plus<double>());
    for(int i=0;i<numPatterns;i++)
        sum[i]/=population_num;
    return sum;
}

se::dd1 se::cal_mean(dd3& all_sol){
    dd1 sum(numPatterns,0.0);
    int population_num=all_sol.size();
    for(int i=0;i<region_num;i++){
        for(int j=0;j<sample_num;j++)
            transform(sum.begin(),sum.end(),all_sol[i][j].begin(),sum.begin(),plus<double>());
    }
    for(int i=0;i<numPatterns;i++)
        sum[i]/=(region_num*sample_num);
    return sum;
}
void se::cal_ev(){

    double Tj,Vj,Mj;
    for(int i=0;i<searcher_num;i++)
        searcher_fit[i]=evaluation(searcher[i]);
    //計算M_j
    double sample_fit_sum=0.0;
    for(int i=0;i<region_num;i++){
        //存取以往最好的rbj---------------------------------------------------------------------------------------------(要修)//
        double tmp_best=region_best_fit[i];
        int best_pos=-1;

        for(int j=0;j<sample_num;j++){
            if(not_evaluate)
                sample_fit[i][j]=evaluation(sample[i][j]);
            // sample_fit_sum+=1/sample_fit[i][j];//越小越好 所以取倒數
            sample_fit_sum+=sample_fit[i][j];
            if(sample_fit[i][j]<tmp_best){
                tmp_best=sample_fit[i][j];
                best_pos=j;
            }
        }
        if (best_pos >= 0) {
            region_best_fit[i] = tmp_best;
            region_best[i]=sample[i][best_pos];
        }
    }
    //計算V_j
    dd2 Vij(searcher_num,dd1(region_num,0.0));

     for (int i = 0; i < searcher_num; i++) {
        for (int j = 0; j < region_num; j++) {
            for (int k = 0; k < sample_num; k++) {
                sampleV_fit[i][j][k] = evaluation(sampleV[i][j][k]);
                // Vij[i][j] += (1/sampleV_fit[i][j][k]);//越小越好 所以取倒數
                Vij[i][j] += (sampleV_fit[i][j][k]);
            }
            Vij[i][j] /= sample_num; 
        }
    }

     for (int i = 0; i < searcher_num; i++) 
        for (int j = 0; j < region_num; j++) 
            EV[i][j] = (tb[j]/ta[j])*(1.0/Vij[i][j])*(1.0/(region_best_fit[j]/sample_fit_sum));
        
    
    choose_sample();
}
//更新sample
void se::choose_sample(){
    dd2 best_sampleV_fit(region_num,dd1(sample_num,inf));
    dd3 best_sampleV(region_num,dd2(sample_num,dd1(numPatterns)));

     for (int i = 0; i < searcher_num ; i++){
        for(int j = 0; j < region_num; j++) {
            for (int k = 0; k < sample_num; k++){
                /*
                if (sampleV_fit[i][j][k] < sample_fit[j][k]) {
                    sample[j][k]= sampleV[i][j][k];
                    sample_fit[j][k] = sampleV_fit[i][j][k];
                }
                */
            //    cout<<"Searcher "<<i<<" Region : "<<j<<" Sample:"<<k<<" fit="<<sampleV_fit[i][j][k]<<endl;
               if(best_sampleV_fit[j][k]>sampleV_fit[i][j][k]){
                   best_sampleV_fit[j][k]=sampleV_fit[i][j][k];
                   best_sampleV[j][k]=sampleV[i][j][k];
               }
            }
        }
        // cout<<endl;
    }

    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++)
        {
            if((double)rand()/RAND_MAX>0.5){
                sample_fit[i][j]=best_sampleV_fit[i][j];
                sample[i][j]=best_sampleV[i][j];
            }
            else{
                int rand_searcher=rand()%searcher_num;
                sample_fit[i][j]=sampleV_fit[rand_searcher][i][j];
                sample[i][j]=sampleV[rand_searcher][i][j];
            }
            // cout<<"Region: "<<i<<" Sample "<<j<<" fit = "<<sample_fit[i][j]<<endl;
        }
    
    
    for(int i=0;i<searcher_num;i++)
        sampleV[i]=sample;
    
}
//兩種版本測試:1.無競技選擇(直接挑最好的),2.有競技選擇
//目前使用2
void se::select_player(){
    for (int i = 0; i < region_num; i++)
        tb[i]++;

    for (int i = 0; i < searcher_num; i++) {
        //使用tournament

        int first_select = rand() % region_num;
        double first_ev = EV[i][first_select];
        for(int j=0;j<tour_size;j++){
            int second_select= rand() % region_num;
            if (EV[i][second_select] >  first_ev) {
                first_select = second_select;
                first_ev = EV[i][second_select];
            }
        }
        
        //不使用tournament
        /*
        int first_select=0;
        for(int j=1;j<region_num;j++)
            if(EV[i][first_select]<EV[i][j])
                first_select=j;
        */
        //cout<<"searcher "<<i<<" select "<<first_select<<" EV="<<EV[i][first_select]<<endl;//------------------------------------//
        int best_sam_pos=0;
        bool flag=false;
        for (int j = 0; j < sample_num; j++){
            // if (sample_fit[first_select][j] < searcher_fit[i]) {
            //     // searcher[i] = sample[first_select][j];
            //     // searcher_fit[i] = sample_fit[first_select][j];
            //     best_sam_pos=j;
            //     flag=true;
            // }
            searcher[i]=region_best[first_select];
        }
        // if(flag){
        //     searcher[i] = sample[first_select][best_sam_pos];
        //     searcher_fit[i] = sample_fit[first_select][best_sam_pos];
        // }
        // update region_it[i] and region_hl[i];
        ta[first_select]++;
        tb[first_select] = 1;
    }
}
se::dd1 se::run(){
    double avg_obj_value;
    double E1;
    gnu_eva=dd1(max_evaluate,0.0);
    vector<double>iter_obj_avg(numIter,0.0);
    for(int i=0;i<numRuns;i++){
        init();
        arrange();
        not_evaluate=1;
        avg_obj_value=0.0;
        for(int j=0;j<numIter;j++){
            E1=2.0*(1.0-((double)j/(double)numIter));
            HHO_transit(E1);
            cal_ev();
            select_player();
            marketing_survey();
            
            not_evaluate=0;
            iter_obj_avg[j]+=best_obj_val;
        }
        
    }
    ofstream output("gnuplot/result/SEHHO.txt");
    for(int i=0;i<max_evaluate;i++)
        if(i%500==0){
            output<<i<<" ";
            output<<fixed<<setprecision(20)<<gnu_eva[i]/numRuns<<endl;
        }
    output.close();
    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(20)<<iter_obj_avg[i]/numRuns<<endl;
    cout<<"Evaluation Number:"<<evaluate_num<<endl;
    cout<<"Best objective value:"<<evaluation(best_sol)<<endl;
    return best_sol;

}

void se::marketing_survey()
{
    for (int i = 0; i < region_num; i++)

        if (tb[i] > 1)
            tb[i] = 1.0;

    int j = -1;
    for (int i = 0; i < searcher_num; i++) 
        if (searcher_fit[i] < best_obj_val) {
            best_obj_val = searcher_fit[i];
            j = i;
        }
    
    if (j >= 0)
        best_sol = searcher[j];
}

double se::evaluation(dd1& sol){
    double sum=0.0;
    for(int i=0;i<numPatterns-1;i++)
        sum+=100.0*pow(sol[i+1]-pow(sol[i],2.0),2.0)+pow(sol[i]-1,2.0);
    evaluate_num++;
    gnu_eva[evaluate_num]+=best_obj_val;
    return sum;
}
#endif
       