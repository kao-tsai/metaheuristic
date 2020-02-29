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
using namespace std;

class se{
	public:
		typedef vector<int> d1;
        typedef vector<vector<int>> d2;
        typedef vector<vector<vector<int>>> d3;
        typedef vector<double> dd1;
        typedef vector<vector<double>> dd2;
        typedef vector<vector<vector<double>>> dd3;

		se(int,int,int,int,int,int,int);
		d1 run();
	private:
		void init();
		int evaluation(d1&);
        void ga_crossover();
        void ga_mut();
        void arrange();
        void cal_ev();
        void choose_sample();
        void select_player();
        void marketing_survey();

private:
	int numRuns;
	int numIter;
	int numPatterns;
    
    d2 searcher;
    d3 sample;
    vector<d3> sampleV;
    d2 region_bit;
    dd1 ta;
    dd1 tb;
    d1 searcher_belong;
    dd2 EV;
    dd1 region_best_fit;
    dd2 sample_fit;
    dd3 sampleV_fit;
    dd1 searcher_fit;
    d1 best_sol;
    int best_obj_val;
    int clip_bit_num;

    int region_num;
    int searcher_num;
    int sample_num;
    int tour_size;
};

se::se(int xNumRuns,
int xNumIter,
int xNumPatterns,
int xRegion,
int in_searcher_num,
int in_sample_num,
int in_tournament_size
)
{
    srand(time(0));

    numRuns=xNumRuns;
	numIter=xNumIter;
	numPatterns=xNumPatterns;
    // crossover_pro=0.6;
    // mutation_pro=0.1;
    region_num=xRegion;
    searcher_num=in_searcher_num;
    sample_num=in_sample_num;
    tour_size=in_tournament_size;
}


void se::init(){
    best_obj_val=0;
    best_sol.assign(numPatterns,0);
    tb.assign(region_num,0.0);
    ta.assign(region_num,1.0);
    searcher_belong.assign(searcher_num,0);
    searcher.assign(searcher_num,d1(numPatterns,0));
    sample.assign(region_num,d2(sample_num,d1(numPatterns,0)));
    sampleV.assign(searcher_num, d3(region_num, d2(sample_num*2, d1(numPatterns, 0))));
    sample_fit.assign(region_num, dd1(sample_num, 0.0));
    sampleV_fit.assign(searcher_num,dd2(region_num,dd1(sample_num,0.0)));
    searcher_fit.assign(searcher_num,0.0);
    region_best_fit.assign(region_num,0);
	for (int i = 0; i < searcher_num; i++)
        for(int j=0;j<numPatterns;j++)
            searcher[i][j]=rand()%2;
}

void se::arrange(){
    clip_bit_num=log2(region_num);
    region_bit.assign(region_num,d1(clip_bit_num,0));
    for(int i=0;i<region_num;i++){
        int n=clip_bit_num-1;
        int tmp=i;
        while(n>=0){
            region_bit[i][n]=tmp%2;
            tmp/=2;
            n--;
        }
    }
   
    for(int i=0;i<searcher_num;i++){
        searcher_belong[i]=i%region_num;
        for(int j=0;j<clip_bit_num;j++)
            searcher[i][j]=region_bit[searcher_belong[i]][j];
    }
    for(int i=0;i<region_num;i++)
        for(int j=0;j<sample_num;j++)
            for(int k=0;k<numPatterns;k++)
                if(k<clip_bit_num)
                    sample[i][j][k]=region_bit[i][k];
                else
                    sample[i][j][k]=rand()%2;
    for(int i=0;i<searcher_num;i++){
        tb[searcher_belong[i]]++;
        ta[searcher_belong[i]]=1.0;
    }
                
    EV.assign(searcher_num,dd1(region_num,0.0));
    //
    //
    //
}


void se::ga_crossover(){
    int crossover_point,s_pos;
    for(int i=0;i<searcher_num;i++)
        for(int j=0;j<region_num;j++)
            for(int k=0;k<sample_num;k++){
                crossover_point=rand()%(numPatterns-1)+1;
                for (int l = 0; l < clip_bit_num; l++) {
                    sampleV[i][j][k*2][l] = region_bit[j][l];
                    sampleV[i][j][k*2+1][l] = region_bit[j][l];
                }
                for (int l = clip_bit_num; l < numPatterns; l++) {
                    if (l < crossover_point) {
                        sampleV[i][j][k*2][l] = searcher[i][l];
                        sampleV[i][j][k*2+1][l] = sample[j][k][l];
                    }
                    else {
                        sampleV[i][j][k*2][l] = sample[j][k][l];
                        sampleV[i][j][k*2+1][l] = searcher[i][l];
                    }
                }
            }
}
void se::ga_mut(){
    for (int i = 0; i < searcher_num; i++) {
        for (int j = 0; j < region_num; j++) {
            for (int k = 0; k < 2*sample_num; k++) {
                int m = rand() % numPatterns; 
                if (m >= clip_bit_num)
                    sampleV[i][j][k][m] = !sampleV[i][j][k][m];
            }
        }
    }
}
void se::cal_ev(){

    double Tj,Vj,Mj;

    double sample_fit_sum=0.0;
    for(int i=0;i<region_num;i++){
        double tmp_best=0;
        int best_pos=0;
        for(int j=0;j<sample_num;j++){
            sample_fit[i][j]=evaluation(sample[i][j]);
            sample_fit_sum+=sample_fit[i][j];
            if(sample_fit[i][j]>tmp_best){
                tmp_best=sample_fit[i][j];
                best_pos=j;
            }
            region_best_fit[j]=tmp_best;
            
        }
    }
    // for(int i=0;i<region_num;i++)
    //     Mj=region_best_fit[i]/sample_fit_sum;
    dd2 Vij(searcher_num,dd1(region_num,0.0));

     for (int i = 0; i < searcher_num; i++) {
        for (int j = 0; j < region_num; j++) {
            for (int k = 0; k < sample_num; k++) {
                sampleV_fit[i][j][k*2] = evaluation(sampleV[i][j][k*2]);
                sampleV_fit[i][j][k*2+1] = evaluation(sampleV[i][j][k*2+1]);
                Vij[i][j] += sampleV_fit[i][j][k*2] + sampleV_fit[i][j][k*2+1];
            }
            Vij[i][j] /= 2*sample_num; 
        }
    }

     for (int i = 0; i < searcher_num; i++) 
        for (int j = 0; j < region_num; j++) 
            EV[i][j] = (tb[j]/ta[j])*Vij[i][j]*(region_best_fit[j]/sample_fit_sum);
        
    for(int i=0;i<searcher_num;i++)
        searcher_fit[i]=evaluation(searcher[i]);
    choose_sample();
}
void se::choose_sample(){
     for (int i = 0; i < searcher_num ; i++) {
        for(int j = 0; j < region_num; j++) {
            for (int k = 0; k < sample_num; k++) {
                if (sampleV_fit[i][j][k*2] > sample_fit[j][k]) {
                    for (int l = clip_bit_num; l < numPatterns; l++)
                        sample[j][k][l] = sampleV[i][j][k*2][l];
                    sample_fit[j][k] = sampleV_fit[i][j][k*2];
                }
                if (sampleV_fit[i][j][k*2+1] > sample_fit[j][k]) {
                    for (int l = clip_bit_num; l < numPatterns; l++)
                        sample[j][k][l] = sampleV[i][j][k*2+1][l];
                    sample_fit[j][k] = sampleV_fit[i][j][k*2+1];
                }
            }
        }
    }
}
void se::select_player(){
    for (int i = 0; i < region_num; i++)
        tb[i]++;

    for (int i = 0; i < searcher_num; i++) {

        int first_select = rand() % region_num;
        double first_ev = EV[i][first_select];
        
      
        int second_select= rand() % region_num;
        if (EV[i][second_select] >  first_select) {
            first_select = second_select;
            first_ev = EV[i][first_select];
        }
        


        for (int j = 0; j < sample_num; j++) {
            if (sample_fit[first_select][j] > searcher_fit[i]) {
                searcher[i] = sample[first_select][j];
                searcher_fit[i] = sample_fit[first_select][j];
            }
        }

        // update region_it[i] and region_hl[i];
        ta[first_select]++;
        tb[first_select] = 1;
    }
}
se::d1 se::run(){

    dd1 iter_obj_avg(numIter,0.0);
    for(int i=0;i<numRuns;i++){
        //cout<<"hello:"<<count<<endl;

        init();
        
        arrange();
        for(int j=0;j<numIter;j++){
            
            ga_crossover();
            ga_mut();
            cal_ev();
            select_player();
            
            marketing_survey();
            iter_obj_avg[j]+=best_obj_val;
        }
        
    }

    for(int i=0;i<numIter;i++)
        cout<<fixed<<setprecision(3)<<iter_obj_avg[i]/numRuns<<endl;

    return best_sol;

}

void se::marketing_survey()
{
    for (int i = 0; i < region_num; i++)

        if (tb[i] > 1)
            tb[i] = 1.0;

    int j = -1;
    for (int i = 0; i < searcher_num; i++) {
        if (searcher_fit[i] > best_obj_val) {
            best_obj_val = searcher_fit[i];
            j = i;
        }
    }
    if (j >= 0)
        best_sol = searcher[j];
}

int se::evaluation(d1& sol){
    int  count=0;
    for(int i=0;i<numPatterns;i++)
	    count+=sol[i];
    return count;
}
#endif
       