#ifndef __MOPSO_H_INCLUDED__
#define __MOPSO_H_INCLUDED__
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include "../test_problem/dataset.cpp"
#include "../test_problem/function_init.cpp"

using namespace std;
# define EPS 1.0e-14
class mopso{
    public:
		typedef vector<double> solution;
        typedef vector<solution> population;
		mopso(int,int,int,double,int,string);
		population run();
	private:
		void init(population&);
        population fitness(population&);
        vector<population> fast_non_dominated(vector<population>&,vector<vector<double>>);
        void obj_sorting(solution&,vector<int>&);
        void update_pbest(population&);
        void update_cube(population &);
         int cal_cube_id(solution&,solution&,solution&);
        void init_min_max();
        population get_non_set(population&,population&);
        bool dominate(solution&,solution&);
        void create_cube();
        void delete_par();
        int select_lead();
        void update_wheel();
        void update_velocity(population&,population&);

private:
	int numRuns;
	int numIter;
	population obj_val;
    int population_num;
    string problem_func;
    int dimension; //每一個解的維度
    int obj_num; //目標數量
    double ndiv;
    int mem;
    solution upperbound;
    solution lowerbound;
    vector<int> cube_id;
    //vector<unsigned int>count_id;
    population hypercube;
    population hyper_fit;
    solution hyper_max;
    solution hyper_min;
    population pbest;
    population pbest_fit;
    vector<vector<int>> wheel;
};
mopso::mopso(int xNumRuns,
int xNumIter,
int xPopulationNum,
double x_ndiv,
int x_mem,
string xproblem_func
)
{
    srand(time(0));
    numRuns=xNumRuns;
	numIter=xNumIter;
    population_num=xPopulationNum;
    problem_func=xproblem_func;
    func_init(problem_func);
    dimension=idimension;
    obj_num=iobj_num;	
    upperbound=ub;
    lowerbound=lb;
    ndiv=x_ndiv;
    mem=x_mem;
   
}

void mopso::init(population& all_sol){
    for(int i=0;i<population_num;i++)
        for(int j=0;j<dimension;j++)
            all_sol[i][j]=(upperbound[j]-lowerbound[j])*((double)rand()/RAND_MAX)+lowerbound[j];
    
    pbest=all_sol;
    pbest_fit=fitness(pbest);

    hyper_max=solution(obj_num);
    hyper_min=solution(obj_num);

    hyper_fit=pbest_fit;
    hypercube=get_non_set(pbest,hyper_fit);
    
    create_cube();
    
    update_wheel();
    delete_par();
}

void mopso::update_pbest(population &new_sol){
    population new_fit;
    new_fit=fitness(new_sol);
    for(int i=0;i<population_num;i++)
        if(dominate(new_fit[i],pbest_fit[i])){
            pbest[i]=new_sol[i];
            pbest_fit[i]=new_fit[i];
        }
}


//-------------------------------找hypercube內各目標的最大和最小值--------------------------------------//
void mopso::init_min_max(){
    for(int i=0;i<obj_num;i++){
        hyper_max[i]=hyper_fit[0][i];
        hyper_min[i]=hyper_fit[0][i];
    }
    for(int i=0;i<obj_num;i++)
        for(int j=1;j<population_num;j++)
            if(hyper_max[i]<hyper_fit[j][i])
                hyper_max[i]=hyper_fit[j][i];
            else if(hyper_min[i]>hyper_fit[j][i])
                hyper_min[i]=hyper_fit[j][i];
}
//---------------------------------------------------------------------------------------------------------------------//
//----------------------------------創建新的HyperCube並分配Cube的ID-----------------------------------//
void mopso::create_cube(){
    init_min_max();
    
    solution dc(obj_num);
    solution start(obj_num);
   
    for(int i=0;i<obj_num;i++){
        //兩種寫法：ndiv減一或不減
         
        dc[i]=(hyper_max[i]-hyper_min[i])/(double)(ndiv-1);
        
        start[i]=hyper_min[i]-(dc[i]/2.0);

    }
    
    int hypercube_size=hypercube.size();
    //count_id=vector<unsigned int>(max_cube,0);
    cube_id=vector< int>(hypercube_size);
    for(int i=0;i<hypercube_size;i++){
        cube_id[i]=cal_cube_id(hyper_fit[i],dc,start);
        //count_id[cube_id[i]]++;
    }
}
//--------------------------------------------------------------------------------------------------------------------//
 int mopso::cal_cube_id(solution& fit,solution& dc, solution& start){
     int id=0,id_dim;
    int j=obj_num-1;
    for(int i=0;i<obj_num;i++,j--){
        id_dim=floor((fit[i]-start[i])/dc[i]);
        id+=floor(id_dim*pow(ndiv,j));
    }
    return id;
}
bool mopso::dominate(solution& new_fit,solution& old_fit){
    if(new_fit[0]<=old_fit[0]&&new_fit[1]<=old_fit[1]){
                if(new_fit[0]==old_fit[0]&&new_fit[1]==old_fit[1])
                    return false;
                else{
                    return true;
                }
    }

}
//------------------------------------------取得Non-Dominated Set---------------------------------------//
mopso::population mopso::get_non_set(population& all_sol,population& fit){
    // all_sol=population(5,solution(1));
    // all_sol[0][0]=1;
    // all_sol[1][0]=-1;
    // all_sol[2][0]=2;
    // all_sol[3][0]=0.5;
    // all_sol[4][0]=3;
    // fit=fitness(all_sol);
    int sol_size=all_sol.size();
    //+sol_size
    vector<bool>check_dominate(sol_size,true);
    //+
    for(int i=0;i<sol_size-1;i++)
        for(int j=i+1;j<sol_size;j++){
            if(dominate(fit[i],fit[j]))
                check_dominate[j]=false;
            else if(dominate(fit[j],fit[i]))
                check_dominate[i]=false;
        }
    //記住山--------------//
    // for(int i=0;i<5;i++){
    //     cout<<i<<":"<<check_dominate[i]<<endl;
    // }
    //------------------------//
    
    population new_non_dom,new_fit;
    for(int i=0;i<sol_size;i++)
        if(check_dominate[i]==true){
            new_non_dom.push_back(all_sol[i]);
            new_fit.push_back(fit[i]);
        }
    fit=new_fit;
    return  new_non_dom;
} 
//-----------------------------------------------------------------------------------------------------------------//
//----------------------------------------Update Wheel Function------------------------------------------//
void mopso::update_wheel(){
    
    vector< int> add_id(1);
    wheel.clear();
    bool flag;
    for(int i=0;i<hypercube.size();i++){
        flag=true;
        for(int j=0;j<wheel.size();j++){
            if(wheel[j][0]==cube_id[i]){
                wheel[j].push_back(i);
                flag=false;
                break;
            }   
        }
        if(flag){
            add_id[0]=cube_id[i];
            wheel.push_back(add_id);
            wheel[wheel.size()-1].push_back(i);
        }
    }
    //記住山//
    // for(int i=0;i<wheel.size();i++)
    //     cout<<"wheel "<<i<<" Size:"<<wheel[i].size()<<endl;

}
//-----------------------------------------------------------------------------------------------------------------//
//------------------------------------------------更新cube----------------------------------------------------//
void mopso::update_cube(population& all_sol){
    population new_nondominated;
    new_nondominated=hypercube;
    new_nondominated.insert(new_nondominated.end(),all_sol.begin(),all_sol.end());
    population tmp_fit;
    tmp_fit=fitness(new_nondominated);
    hypercube=get_non_set(new_nondominated,tmp_fit);
    hyper_fit=tmp_fit;
    
    create_cube();
    
    //------Update Wheel Function------//
    update_wheel();
   
    //當hypercube大小超過上限，從擁擠的hypercube刪除particle//
    delete_par();

    // cout<<"deep"<<endl;
}
//----------------------------------------------------------------------------------------------------------------//

//-------------------------------------刪除HyperCube多出來的粒子------------------------------------//
void mopso::delete_par(){
    
    int tmp_max,rnd1,rnd2;
    vector<int> same_size;
    // int count=0;
    // cout<<"HyperCube size:"<<hypercube.size()<<endl;
    while(hypercube.size()>mem){

        //-----選出有最大值的hypercube id----//
        
        same_size.clear();
        tmp_max=wheel[0].size()-1;
        same_size.push_back(0);
        for(int i=1;i<wheel.size();i++){
            if(tmp_max<wheel[i].size()-1){
                tmp_max=wheel[i].size()-1;
                same_size.clear();
                same_size.push_back(i);
            }
            else if(tmp_max==wheel[i].size()-1)
                same_size.push_back(i);
        }
        //-------------------------------------------------//
        //---rnd1會選出一個相同數目的wheel---//
        rnd1=same_size[rand()%same_size.size()];
        rnd2=rand()%(wheel[rnd1].size()-1)+1;
        
        //刪除 cube_id , hypercube , hyper_fit , 最後再刪除wheel//
        // cout<<"hypercube Size:"<<hypercube.size()<<" cube ID Size:"<<cube_id.size()<<endl;
        // cout<<"rnd1="<<rnd1<<", rnd2="<<rnd2<<endl;
        // for(int i=0;i<wheel.size();i++)
        //     cout<<"wheel "<<i<<" Size:"<<wheel[i].size()<<endl;
        // cout<<wheel[rnd1][rnd2]<<endl;
        cube_id.erase(cube_id.begin()+wheel[rnd1][rnd2]);
        // cout<<"deep1-----------------------------------------------------------------"<<endl;
        hypercube.erase(hypercube.begin()+wheel[rnd1][rnd2]);
        //判斷是否需要刪除hyper_fit
        //hyper_fit.erase(hyper_fit.begin()+wheel[rnd1][rnd2]);
        if(wheel[rnd1].size()==2)
            wheel.erase(wheel.begin()+rnd1);
        else
            wheel[rnd1].erase(wheel[rnd1].begin()+rnd2);
        
        //-----------------------------------------------------------------------//
        // count++;
        // cout<<"deep2:"<<count<<endl;
        update_wheel();
    }
    
}
//---------------------------------------------------------------------------------------------------------------//
int mopso::select_lead(){
    //利用輪盤選擇gbest
    solution prob(wheel.size());
    double hypersum=0;
    for(int i=0;i<wheel.size();i++){
        prob[i]=10.0/(pow(wheel[i].size()-1,3.0));
        // prob[i]=10.0/(wheel[i].size()-1);
        hypersum+=prob[i];
    }
    double rnd1=((double)rand()/RAND_MAX)*hypersum;
    double sum=0;
    int choose;
    for(int i=0;i<wheel.size();i++){
        sum+=prob[i];
        if(rnd1<sum){
            choose=i;
            break;
        }
    }
    int rnd2=rand()%(wheel[choose].size()-1)+1;
    return wheel[choose][rnd2];
}

void mopso::update_velocity(population &all_sol,population &vel){
    double W=0.4,c1=1,c2=2;
    int gbest_pos;
    for(int i=0;i<population_num;i++){
        gbest_pos=select_lead();
        for(int j=0;j<dimension;j++){
            vel[i][j]=W*vel[i][j]+c1*((double)rand()/RAND_MAX)*(pbest[i][j]
                            -all_sol[i][j])+c2*((double)rand()/RAND_MAX)*(hypercube[gbest_pos][j]-all_sol[i][j]);
            all_sol[i][j]+=vel[i][j];
        }
    }
     for(int i=0;i<population_num;i++)
        for(int j=0;j<dimension;j++)
            if(all_sol[i][j]>upperbound[j])
                all_sol[i][j]=upperbound[j];
            else if(all_sol[i][j]<lowerbound[j])
                all_sol[i][j]=lowerbound[j];
        
}


mopso::population mopso::run(){
   
    population currentSol(population_num,solution(dimension));
    population vel;
    for(int i=0;i<numRuns;i++){
        init(currentSol);
        
        vel=population(population_num,solution(dimension,0));
        for(int j=0;j<numIter;j++){
          
            update_velocity(currentSol,vel);
        
            update_pbest(currentSol);
            
            update_cube(currentSol);
            
        }
    }
    hyper_fit=fitness(hypercube);
    return hyper_fit;

}
mopso::population mopso::fitness(population& all_sol) {
    if(problem_func=="SCH")
        return SCH(all_sol);
    else if(problem_func=="FON")
        return FON(all_sol);
    else if(problem_func=="POL")
        return POL(all_sol);
    else if(problem_func=="KUR")
        return KUR(all_sol);
    else if(problem_func=="ZDT1")
        return ZDT1(all_sol);
    else if(problem_func=="ZDT2")
        return ZDT2(all_sol);
    else if(problem_func=="ZDT3")
        return ZDT3(all_sol);
    else if(problem_func=="ZDT4")
        return ZDT4(all_sol);
    else if(problem_func=="ZDT6")
        return ZDT6(all_sol);
    else if(problem_func=="UF1")
        return UF1(all_sol);

}

#endif