#include "test_problem.h"
#define PI  3.1415926535897932384626433832795
using namespace std;
test_problem::test_problem(string func_name,int x_dim){
    if(func_name=="Ackley"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -32;
            ub[i] = 32;
        }
    }
    else if(func_name=="Rastrigin"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -5.12;
            ub[i] = 5.12;
        }
    }
    else if(func_name=="Sphere"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -5.12;
            ub[i] = 5.12;
        }
    }
    else if(func_name=="Rosenbrock"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -5;
            ub[i] = 10;
        }
    }
    else if(func_name=="Michalewicz"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = 0;
            ub[i] = M_PI;
        }
    }
    else if(func_name=="F15"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -20;
            ub[i] = 20;
        }
    }
    else if(func_name=="F16"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -10;
            ub[i] = 10;
        }
    }
    

}

// double test_problem::Ackley(solution1 &sol){
    
//     double sum_of_sq=0.0;
//     double sum_of_cos=0.0;
//     double c=2*PI;
//     for(int i=0;i<idim;i++)
//         sum_of_sq +=pow(sol[i],2.0);
//     for(int i=0;i<idim;i++)
//         sum_of_cos +=cos(c*sol[i]);
    
//     return (-20)*exp((-0.2)*sqrt((1.0/idim)*sum_of_sq))-exp((1.0/idim)*sum_of_cos)+exp(1)+20;
    
// }

double test_problem::Ackley(solution1 &sol){
    double temp1=0,temp2=0;
    for(int d=0;d<idim;d++)
        temp1 += pow(sol[d],2);
    temp1 /= (double)idim;
    temp1 = (double)sqrt(temp1);
    temp1 *= (double)-0.2;

    for(int d=0;d<idim;d++)
        temp2 += (double)cos(2*M_PI*sol[d]);
    temp2 /= idim;

    return (double)-20*exp(temp1)-exp(temp2)+exp(1)+20;
}

double test_problem::Rastrigin(solution1 &sol){
    double sum=0.0;
    for(int i=0;i<idim;i++)
        sum+=pow(sol[i],2.0)-10*cos(2*M_PI*sol[i])+10;
    
    return sum;
}
double test_problem::Sphere(solution1 &sol){
    double sum=0.0;
    for(int i=0;i<idim;i++)
        sum+=pow(sol[i],2.0);
    return sum;
}
double test_problem::Rosenbrock(solution1 &sol){
    double sum=0.0;
    for(int i=0;i<idim-1;i++){
        sum+=100.0*pow((sol[i+1]-pow(sol[i],2.0)),2.0)+pow(sol[i]-1.0,2.0);
    }
    return sum;
}
double test_problem::Michalewicz(solution1 &sol){
    double sum=0.0;
    for(int i=0;i<idim;i++){
        sum+=sin(sol[i])*pow(sin(i*pow(sol[i],2.0)/M_PI),20);
    }
    return -sum;
}
double test_problem::F15(solution1 &sol){
    return 0.0;
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