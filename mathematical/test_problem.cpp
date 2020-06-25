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
    else if(func_name=="Griewank10"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -600;
            ub[i] = 600;
        }
    }else if(func_name=="Schaffer2"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -100;
            ub[i] = 100;
        }
    }else if(func_name=="Schwefel"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -10;
            ub[i] = 10;
        }
    }else if(func_name=="Bohachevsky1"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -50;
            ub[i] = 50;
        }
    }else if(func_name=="Sum_Square"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -10;
            ub[i] = 10;
        }
    }else if(func_name=="Booth"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -10;
            ub[i] = 10;
        }
    }else if(func_name=="Zakharov"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -5;
            ub[i] = 10;
        }
    }else if(func_name=="Three_Hump_Camel"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -5;
            ub[i] = 5;
        }
    }else if(func_name=="De_Jong5"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -65.536;
            ub[i] = 65.536;
        }
    }else if(func_name=="Beale"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -4.5;
            ub[i] = 4.5;
        }
    }else if(func_name=="Powell"){
        idim=x_dim;
        ub.assign(idim, 0);
        lb.assign(idim, 0);
        for(int i=0;i<idim;i++){
            lb[i] = -4;
            ub[i] = 5;
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
double test_problem::Griewank10(solution1& Xi){ // Griewank2 [-600, 600]
    double fit=0,temp1=0,temp2=1;
    for(int d=0; d<idim; ++d){
        temp1 += pow(Xi[d]-100, 2);
        temp2 *= cos((Xi[d]-100)/(sqrt(d+1)));
    }
    fit = (double)temp1/4000 - temp2 + 1;
    return fit;
}
double test_problem::Schaffer2(solution1 &Xi){ // Schaffer2 [-100,100]
    double fit=0;
    fit = pow(Xi[0]*Xi[0]+Xi[1]*Xi[1], 0.25) * (50 * pow(Xi[0]*Xi[0]+Xi[1]*Xi[1], 0.1) + 1);
    return fit;
}
double test_problem::Schwefel(solution1 &Xi){ // Schwefel [-10, 10]
    double fit=0, temp=0;
    for(int d=0; d<idim; d++){
        temp = 0;
        for(int i=0; i<d; i++)
            temp += Xi[d];
        fit += temp * temp;
    }
    return fit;
}
double test_problem::Bohachevsky1(solution1 &Xi){ // Bohachevsky [-50, 50]
    double fit=0;
    fit = Xi[0]*Xi[0] + 2*Xi[1]*Xi[1] - 0.3*cos(3*M_PI*Xi[0]) - 0.4*cos(4*M_PI*Xi[1]) + 0.7;
    return fit;
}
double test_problem::Sum_Square(solution1 &Xi){ // SumSquares [-10, 10]
    double fit=0;
    for(int d=0; d<idim; d++)
        fit += (d+1)*Xi[d]*Xi[d];
    return fit;
}
double test_problem::Booth(solution1 &Xi){ // Booth [-10, 10]
    double fit=0;
    fit = pow(Xi[0]+ 2*Xi[1] -7, 2) + pow(2*Xi[0]+Xi[1]-5, 2);
    return fit;
}
double test_problem::Zakharov(solution1 &Xi){ // Zakharov [-5, 10]
    double fit=0;
    double temp1=0, temp2=0;
    for(int d=0; d<idim; d++)
        fit += pow(Xi[d], 2);

    for(int d=0; d<idim; d++)
        temp1 += 0.5*(d+1)*Xi[d];
    temp1 = pow(temp1, 2);

    for(int d=0; d<idim; d++)
        temp2 += 0.5*(d+1)*Xi[d];
    temp2 = pow(temp2, 4);

    fit += temp1 + temp2;

    return fit;
}
double test_problem::Three_Hump_Camel(solution1 &Xi){ // Three-hump camel [-5,5]
    double fit=0;
    fit = 2*pow(Xi[0], 2) - 1.05*pow(Xi[0], 4) + pow(Xi[0], 6)/6 + Xi[0]*Xi[1] + pow(Xi[1], 2);
    return fit;
}
double test_problem::De_Jong5(solution1 &Xi){ // De jong function N.5 [-65.536,65.536]
    vector<int> a = vector<int>(5);  
    a[0] = -32; a[1] = -16; a[2] = 0; a[3] = 16; a[4] = 32;
    double fit=0;

    fit = 0.002;
    for(int i=0; i<5; i++)
        for(int j=0; j<5; j++)
            fit += 1/( (i+1) + pow(Xi[0]-a[j%5] ,6) + pow(Xi[1]-a[i%5] ,6) );
        
    fit = 1/fit;

    return fit;
}
double test_problem::Beale(solution1 &Xi){ // Beale [-4.5,4.5]
    double fit=0;
    fit = pow(1.5-Xi[0]+Xi[0]*Xi[1], 2) + pow(2.25-Xi[0]+Xi[0]*Xi[1]*Xi[1], 2) + pow(2.625-Xi[0]+Xi[0]*pow(Xi[1], 3), 2);
    return fit;
}
double test_problem::Powell(solution1 &Xi){ // Powell [-4,5]
    double fit=0;
    for(int i=0; i<idim/4; i++)
        fit += pow(Xi[4*i-3]+10*Xi[4*i-2], 2) + 5*pow(Xi[4*i-1]+Xi[4*i], 2) + pow(Xi[4*i-2]-2*Xi[4*i-1], 4) + 10*pow(Xi[4*i-3]-Xi[4*i], 4);
    return fit;
}
