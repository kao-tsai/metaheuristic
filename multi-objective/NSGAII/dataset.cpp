#include<iostream>
#include<vector>
#include<cmath>
using namespace std;
typedef vector<double> solution1;
typedef vector<solution1> population1;

population1 SCH(population1 &all_sol){
    int obj_num=2;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    
    for(int i=0;i<all_sol.size();i++){
        all_sol_obj[i][0]=pow(all_sol[i][0],2.0);
        all_sol_obj[i][1]=pow(all_sol[i][0]-2.0,2.0);
    }
        
    return all_sol_obj;
}

population1 FON(population1 &all_sol){
    int obj_num=2;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        for(int j=0;j<3;j++)
        {
            f1+=pow(all_sol[i][j]-(1/sqrt(3.0)),2.0);
            f2+=pow(all_sol[i][j]+(1/sqrt(3.0)),2.0);
        }
        f1=1-exp(-f1);
        f2=1-exp(-f2);
        all_sol_obj[i][0]=f1;
        all_sol_obj[i][1]=f2;
    }
        
    return all_sol_obj;
}

population1 POL(population1 &all_sol){
    int obj_num=2;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
    double A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
    double B1;
    double B2;
    for(int i=0;i<all_sol.size();i++){
        B1 = 0.5*sin(all_sol[i][0])-2*cos(all_sol[i][0])+sin(all_sol[i][1])-1.5*cos(all_sol[i][1]);
        B2 = 1.5*sin(all_sol[i][0])-cos(all_sol[i][0])+2*sin(all_sol[i][1])-0.5*cos(all_sol[i][1]);
        f1=1+pow(A1-B1,2)+pow(A2-B2,2);
        f2=pow(all_sol[i][0]+3,2)+pow(all_sol[i][1]+1,2);
        all_sol_obj[i][0]=f1;
        all_sol_obj[i][1]=f2;
    }
    return all_sol_obj;
}
population1 KUR(population1 &all_sol){
    int obj_num=2;
    int dimension=3;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
  
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        for(int j=0;j<dimension-1;j++){
            f1 += -10*exp(-0.2*sqrt(pow(all_sol[i][j],2)+pow(all_sol[i][j+1],2)));
        }

        for(int j=0;j<dimension;j++){
            f2 += pow(abs(all_sol[i][j]),0.8) + 5*sin(pow(all_sol[i][j],3));
        }
        all_sol_obj[i][0] = f1;
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}
population1 ZDT1(population1 &all_sol){
    int obj_num=2;
    int dimension=30;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1 = all_sol[i][0];
        all_sol_obj[i][0] = f1;

        for(int j=1;j<dimension;j++){
            g += all_sol[i][j];
        }
        g = 1 + 9 * g / (dimension-1);
        
        f2 = g*(1-sqrt(all_sol[i][0]/g));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}
population1 ZDT2(population1 &all_sol){
    int obj_num=2;
    int dimension=30;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1 = all_sol[i][0];
        all_sol_obj[i][0] = f1;

        for(int j=1;j<dimension;j++){
            g += all_sol[i][j];
        }
        g = 1 + 9 * g / (dimension-1);
        
        f2 = g*(1-pow(all_sol[i][0]/g,2.0));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}
population1 ZDT3(population1 &all_sol){
    int obj_num=2;
    int dimension=30;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1 = all_sol[i][0];
        all_sol_obj[i][0] = f1;

        for(int j=1;j<dimension;j++){
            g += all_sol[i][j];
        }
        g = 1 + 9 * g / (dimension-1);
        
        f2 = g*(1 - sqrt(all_sol[i][0]/g) - (all_sol[i][0]/g)*sin(10*M_PI*all_sol[i][0]));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}
population1 ZDT4(population1 &all_sol){
    int obj_num=2;
    int dimension=10;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        all_sol_obj[i][0] = all_sol[i][0];

        for(int j=1;j<dimension;j++){
            g += pow(all_sol[i][j],2)-10*cos(4*M_PI*all_sol[i][j]);
        }
        g = 1+10*(dimension-1)+g;
        
        f2 = g*(1-sqrt(all_sol[i][0]/g));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}

population1 ZDT6(population1 &all_sol){
    int obj_num=2;
    int dimension=10;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1=1 - exp(-4*all_sol[i][0])*pow(sin(6*M_PI*all_sol[i][0]),6);
        all_sol_obj[i][0]=f1;

       for(int j=1;j<dimension;j++){
            g += all_sol[i][j];
        }
        g = 1+9*sqrt(sqrt(g/(dimension-1)));
        
        f2 = g*(1-pow(f1/g,2));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}
population1 UF1(population1 &all_sol){
    int obj_num=2;
    int dimension=10;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double f1;
    double f2;
    double g;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1=1 - exp(-4*all_sol[i][0])*pow(sin(6*M_PI*all_sol[i][0]),6);
        all_sol_obj[i][0]=f1;

       for(int j=1;j<dimension;j++){
            g += all_sol[i][j];
        }
        g = 1+9*sqrt(sqrt(g/(dimension-1)));
        
        f2 = g*(1-pow(f1/g,2));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}