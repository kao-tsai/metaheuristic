#ifndef __LIB_H__
#define __LIB_H__

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
using namespace std;

typedef vector<double> sol;
typedef vector<sol> pop;

class LIB{ 
    
public:
    int dim;
    int pop_num;
    double upper_bound;
    double lower_bound;

    LIB(int dim, int pop_num, double upper_bound, double lower_bound);
    double F1(const sol &Xi); // Sphere
    double F2(const sol &Xi);
    double F3(const sol &Xi);
    double F4(const sol &Xi);
    double F5(const sol &Xi); // Rosenbrock
    double F6(const sol &Xi);
    double F7(const sol &Xi);
    double F8(const sol &Xi);
    double F9(const sol &Xi); // Rastrigin
    double F10(const sol &Xi); // Ackley
    double F11(const sol &Xi);
    double F12(const sol &Xi);
    double F13(const sol &Xi);
    double F14(const sol &Xi);
    double F15(const sol &Xi);
    double F16(const sol &Xi);
    double F17(const sol &Xi);
    double F18(const sol &Xi);
    double F19(const sol &Xi);
    double F20(const sol &Xi);
    double F21(const sol &Xi);
    double F22(const sol &Xi);
    double F23(const sol &Xi);
    double F24(const sol &Xi);
    double penalized_function(const int Xi, const int a, const int k, const int m);

    void Print_sol(const sol &sol);
    void Print_pop(const pop &pop);

    int Partition(sol &fit, pop &X, int left, int right); // 同時排序 (index,value)
    void QuickSort(sol &fit, pop &X, int left, int right);

    pop Init_pop();
    sol Check_bound(const sol &sol);

};

LIB::LIB(int dim, int pop_num, double upper_bound, double lower_bound):
        dim(dim), pop_num(pop_num), upper_bound(upper_bound), lower_bound(lower_bound){
    srand(time(NULL));
}

double LIB::F1(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d)
        fit += pow(Xi[d],2); 

    return fit;
}
double LIB::F2(const sol &Xi){
    double fit=0,temp=1;
    for(int d=0; d<dim; ++d)
        fit += abs(Xi[d]);
    for(int d=0; d<dim; ++d)
        temp *= abs(Xi[d]);
    fit += temp;

    return fit;
}
double LIB::F3(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d){
        double sum=0;
        for(int i=0; i<d; ++i)
            sum += Xi[i];
        
        fit += pow(sum,2);  
    }
    return fit;
}
double LIB::F4(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d)
        if(fit<Xi[d])
            fit = Xi[d];
    return fit;
}
double LIB::F5(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim-1; ++d)
        fit += 100*pow(Xi[d+1] + pow(Xi[d],2), 2) + pow(Xi[d]-1, 2);
    return fit;
}
double LIB::F6(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d)
        fit += pow(Xi[d]+0.5, 2);
    return fit;
}
double LIB::F7(const sol &Xi){
    double fit=0,rand_num=0;
    rand_num = (double)(rand()%(RAND_MAX-1))/RAND_MAX;
    for(int d=0; d<dim; ++d)
        fit += (d+1)*pow(Xi[d], 4) + rand_num;
    return fit;
}
double LIB::F8(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d)
        fit += -Xi[d] * sin(sqrt(abs(Xi[d])));
    return fit;
}
double LIB::F9(const sol &Xi){ // Rastrigin
    double fit=0;
    for(int d=0; d<dim; ++d)
        fit += pow(Xi[d], 2) - 10*cos(2*M_PI*Xi[d]) + 10;
    return fit;
}
double LIB::F10(const sol &Xi){ // Ackley
    double temp1=0,temp2=0;
    for(int d=0;d<dim;d++)
        temp1 += pow(Xi[d],2);
    temp1 /= (double)dim;
    temp1 = (double)sqrt(temp1);
    temp1 *= (double)-0.2;

    for(int d=0;d<dim;d++)
        temp2 += (double)cos(2*M_PI*Xi[d]);
    temp2 /= dim;

    return (double)-20*exp(temp1)-exp(temp2)+exp(1)+20;
}
double LIB::F11(const sol &Xi){
    double fit=0,temp1=0,temp2=1;
    for(int d=0; d<dim; ++d){
        temp1 += Xi[d];
        temp2 *= cos(Xi[d]/(sqrt(d+1)));
    }
    fit = temp1/4000 - temp2 + 1;
    return fit;
}
double LIB::F12(const sol &Xi){ 
    double fit=0;
    for(int d=0; d<dim-1; ++d){
        double Yi=1+(double)Xi[d]/4; 
        double Ynext=1+(double)Xi[d+1]/4; 
        fit += pow(Yi-1, 2) * (1+ 10*pow(sin(M_PI*Ynext),2));
    }
    double Y1=1+(double)Xi[0]/4;
    double Yn=1+(double)Xi[dim-1]/4;
    fit += 10*sin( M_PI*Y1 ) + pow(Yn-1,2);
    fit = fit * M_PI / dim;
    for(int d=0; d<dim; ++d)
        fit += penalized_function(Xi[d],10,100,4);
    return fit;
}
double LIB::F13(const sol &Xi){
    double fit=0;
    for(int d=0; d<dim; ++d)
        fit += pow(Xi[d]-1, 2) * (1+ pow(sin(3*M_PI*Xi[d]+1),2));
    
    fit += pow(sin(3*M_PI*Xi[0]),2) + pow(Xi[dim-1],2)*( 1+pow(sin(2*M_PI*Xi[dim-1]),2) );
    fit *= 0.1;
    for(int d=0; d<dim; ++d)
        fit += penalized_function(Xi[d],5,100,4);
    return fit;
}
double LIB::F14(const sol &Xi){
    double fit=0;
    for(int i=0;i<dim;i++)
        fit+=sin(Xi[i])*pow(sin(i*pow(Xi[i],2.0)/M_PI),20);
    
    return -fit;
}
double LIB::F15(const sol &Xi){
    double fit=0;
    for(int i=0;i<dim;i++)
        fit+=sin(Xi[i])*pow(sin(i*pow(Xi[i],2.0)/M_PI),20);
    
    return -fit;
}
double LIB::F16(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F17(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F18(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F19(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F21(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F22(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F23(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::F24(const sol &Xi){
    double fit=0;
    return fit;
}
double LIB::penalized_function(const int Xi, const int a, const int k, const int m){
    double fit=0;
    if(Xi>a)
        return k*pow(Xi-a, m);
    else if(Xi<-a)
        return k*pow(-Xi-a, m);
    else
        return 0;
}
void LIB::Print_sol(const sol &sol){
    for(int d=0; d<dim; ++d)
        cout << sol[d] << " ";
    cout << endl;
}
void LIB::Print_pop(const pop &pop){
    for(int p=0; p<pop_num; ++p){
        for(int d=0; d<dim; ++d)
            cout << pop[p][d] << " ";
        cout << endl; 
    }
}

pop LIB::Init_pop(){

    pop population = pop(pop_num,sol(dim));

    for(int p=0; p<pop_num; ++p){
        sol temp_sol = sol(dim);
        for(int d=0; d<dim; ++d){
            temp_sol[d] = (double) (rand()%RAND_MAX)/RAND_MAX;
            temp_sol[d] *= (double) (upper_bound-lower_bound);
            temp_sol[d] += (double) lower_bound;
        }
        population[p] = temp_sol;
    }

    return population;

}

sol LIB::Check_bound(const sol &solution){
    sol temp_sol = solution;

    double gap;

    // for(int d=0; d<dim; ++d){
    //     if(temp_sol[d]>upper_bound)
    //         temp_sol[d] = upper_bound;
    //     else if(temp_sol[d]<lower_bound)
    //         temp_sol[d] = lower_bound;
    // }

    
    for(int d=0; d<dim; ++d){
        if(temp_sol[d]>upper_bound || temp_sol[d]<lower_bound ){
            temp_sol[d] = (double) rand()/RAND_MAX;
            temp_sol[d] *= (double) (upper_bound-lower_bound);
            temp_sol[d] += (double) lower_bound;
        }
    }

    return temp_sol;
}

int LIB::Partition(sol &fit, pop &X, int left, int right){
    double pivot = fit[right];
    int i = left -1;
    for (int j = left; j < right; j++) {
        if (fit[j] < pivot) {
            i++;
            double temp = fit[i];
            fit[i] = fit[j];
            fit[j] = temp;
            X[i].swap(X[j]);
        }
    }
    i++;
    double temp = fit[i];
    fit[i] = fit[right];
    fit[right] = temp;
    X[i].swap(X[right]);
    return i;
}
void LIB::QuickSort(sol &fit, pop &X, int left, int right){
    if (left < right) {
        int pivot = Partition(fit, X, left, right);
        QuickSort(fit, X, left, pivot - 1);
        QuickSort(fit, X, pivot + 1, right);
    }
}
double test_problem::F16(solution1 &sol){
    double sum1=0.0,sum2=0.0,sum3=0.0;
    for(int i=0;i<idim;i++){
        sum1+=pow(sin(sol[i]),2.0);
        sum2+=pow(sol[i],2.0);
        sum3+=pow(sin(sqrt(fabs(sol[i]))),2.0);
    }
    return (sum1-exp(-sum2))*exp(-sum3);
}
#endif
