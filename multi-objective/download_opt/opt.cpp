#include "opt.h"

opt::opt(int num){
    opt_num = num;
}

void opt::SCH(){
    double min = 0;
    double max = 2;
    int dimension = 1;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){
        for(int j=0;j<dimension;j++){
            value[i][j] = value[i-1][j] + (max-min)/(opt_num-1);
        }
    }

    fstream fpareto;
    string fparetoname = "SCH_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        pareto[i][0] = pow(value[i][0], 2);
        pareto[i][1] = pow(value[i][0] - 2, 2);
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;
    }

    fpareto.close();
}

void opt::FON(){
    double min = -1/pow(3, 0.5);
    double max = 1/pow(3, 0.5);
    int dimension = 3;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){
        for(int j=0;j<dimension;j++){
            value[i][j] = value[i-1][j] + (max-min)/(opt_num-1);
        }
    }

    fstream fpareto;
    string fparetoname = "FON_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        double sum1 = 0;
        double sum2 = 0;
        for(int j=0;j<dimension;j++){
            sum1 += pow(value[i][j]-1/sqrt(3),2);
        }
        sum1 *= -1;
        sum1 = 1-exp(sum1);
        for(int j=0;j<dimension;j++){
            sum2 += pow(value[i][j]+1/sqrt(3),2);
        }
        sum2 *= -1;
        sum2 = 1-exp(sum2);

        pareto[i][0] = sum1;
        pareto[i][1] = sum2;
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;
    }

    fpareto.close();
}
/*
void opt::KUR(){
    double min = -1.2;
    double max = -0.6;
    int dimension = 3;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){
        for(int j=0;j<dimension;j++){
            value[i][j] = value[i-1][j] + (max-min)/(opt_num-1);
        }
    }

    fstream fpareto;
    string fparetoname = "KUR_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    double f1;
    double f2;
    for(int i=0;i<opt_num;i++){
        f1=0;
        f2=0;
        for(int j=0;j<dimension-1;j++){
            f1 += -10*exp(-0.2*sqrt(pow(value[i][j],2)+pow(value[i][j+1],2)));
        }

        for(int j=0;j<dimension;j++){
            f2 += pow(abs(value[i][j]),0.8) + 5*sin(pow(value[i][j],3));
        }
        pareto[i][0] = f1;
        pareto[i][1] = f2;
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;
    }
    fpareto.close();
}*/

void opt::POL(){
    double min = 0;
    double max = 2;
    int dimension = 1;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){
        for(int j=0;j<dimension;j++){
            value[i][j] = value[i-1][j] + (max-min)/(opt_num-1);
        }
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        
    }
}

void opt::ZDT1(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    fstream fpareto;
    string fparetoname = "ZDT1_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        pareto[i][0] = value[i][0];        
        pareto[i][1] = 1-sqrt(value[i][0]);
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::ZDT2(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    fstream fpareto;
    string fparetoname = "ZDT2_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        pareto[i][0] = value[i][0];
        pareto[i][1] = 1-pow(value[i][0],2);
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::ZDT3(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    fstream fpareto;
    string fparetoname = "ZDT3_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        pareto[i][0] = value[i][0];
        pareto[i][1] = 1 - sqrt(value[i][0]) - (value[i][0])*sin(10*M_PI*value[i][0]);
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::ZDT4(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    fstream fpareto;
    string fparetoname = "ZDT4_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }

    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        pareto[i][0] = value[i][0];
        pareto[i][1] = 1-sqrt(value[i][0]);
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;
    }
}

void opt::ZDT6(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    fstream fpareto;
    string fparetoname = "ZDT6_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
       pareto[i][0] = 1 - exp(-4*value[i][0]) * pow(sin(6*M_PI*value[i][0]),6);
       pareto[i][1] = 1 - pow(pareto[i][0],2);
       fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF1(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=2;j<=dimension;j++){
            value[i][j-1] = sin(6.0*PI*value[i][0] + j*PI/dimension);
        }
    }

    fstream fpareto;
    string fparetoname = "UF1_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, sum2, yj2;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++){
            yj1 = value[i][j-1] - sin(6.0*PI*value[i][0] + j*PI/dimension);
            yj1 = yj1 * yj1;
            if(j % 2 == 1) {
                sum1 += yj1;
                count1++;
            } 
        }
        pareto[i][0] = (value[i][0] + 2.0 * sum1 / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++){
            yj2 = value[i][j-1] - sin(6.0*PI*value[i][0] + j*PI/dimension);
            yj2 = yj2 * yj2;
            if(j % 2 == 0) {
                sum2 += yj2;
                count2++;
            } 
        }
        
        pareto[i][1] = (1.0 - sqrt(value[i][0]) + 2.0 * sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF2(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }

    for(int i=0;i<opt_num;i++){
        for(int j=2;j<=dimension;j++){
            if(j % 2 == 1) {
                value[i][j-1] = (0.3*pow(value[i][0],2)*cos(24*PI*value[i][0] + (4*j*PI)/dimension)+0.6*value[i][0])*cos(6*PI*value[i][0] + (j*PI)/dimension);
            }
            if(j % 2 == 0) {
                value[i][j-1] = (0.3*pow(value[i][0],2)*cos(24*PI*value[i][0] + (4*j*PI)/dimension)+0.6*value[i][0])*sin(6*PI*value[i][0] + (j*PI)/dimension);
            }             
        }
    }

    fstream fpareto;
    string fparetoname = "UF2_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, sum2, yj2;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            if(j % 2 == 1) {
                yj1 = value[i][j-1]-0.3*value[i][0]*(value[i][0]*cos(24.0*PI*value[i][0]+4.0*j*PI/dimension)+2.0)*cos(6.0*PI*value[i][0]+j*PI/dimension);
                sum1 += yj1*yj1;
                count1++;
            }
        }
        pareto[i][0] = (value[i][0] + 2.0 * sum1 / (double)count1);
        

        sum2  = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            if(j % 2 == 0) {
                yj2 = value[i][j-1]-0.3*value[i][0]*(value[i][0]*cos(24.0*PI*value[i][0]+4.0*j*PI/dimension)+2.0)*sin(6.0*PI*value[i][0]+j*PI/dimension);
                sum2 += yj2*yj2;
                count2++;
            } 
        }
        pareto[i][1] = (1.0 - sqrt(value[i][0]) + 2.0 * sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF3(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=2;j<=dimension;j++){
            value[i][j-1] = pow(value[i][0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
        }
    }

    fstream fpareto;
    string fparetoname = "UF3_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, prod1, pj1, sum2, yj2, prod2, pj2;
        
        sum1 = 0.0;
        count1 = 0;
        prod1 = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj1 = value[i][j-1]-pow(value[i][0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
            pj1 = cos(20.0*yj1*PI/sqrt(j+0.0));
            if (j % 2 == 1) {
                sum1  += yj1*yj1;
                prod1 *= pj1;
                count1++;
            }
        }
        pareto[i][0] = (value[i][0] + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        prod2  = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj2 = value[i][j-1]-pow(value[i][0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
            pj2 = cos(20.0*yj2*PI/sqrt(j+0.0));
            if (j % 2 == 0) {
                sum2  += yj2*yj2;
                prod2 *= pj2;
                count2++;
            } 
            
        }
        pareto[i][1] = (1.0 - sqrt(value[i][0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF4(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=2;j<=dimension;j++){
            value[i][j-1] = sin(6.0*PI*value[i][0]+j*PI/dimension);
        }
    }

    fstream fpareto;
    string fparetoname = "UF4_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, hj1, sum2, yj2, hj2;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            yj1 = value[i][j-1]-sin(6.0*PI*value[i][0]+j*PI/dimension);
            hj1 = fabs(yj1)/(1.0+exp(2.0*fabs(yj1)));
            if (j % 2 == 1) {
                sum1  += hj1;
                count1++;
            }
        }
        pareto[i][0] = (value[i][0] + 2.0*sum1 / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            yj2 = value[i][j-1]-sin(6.0*PI*value[i][0]+j*PI/dimension);
            hj2 = fabs(yj2)/(1.0+exp(2.0*fabs(yj2)));
            if (j % 2 == 0) {
                sum2  += hj2;
                count2++;
            }
        }
        pareto[i][1] = (1.0 - value[i][0]*value[i][0] + 2.0*sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF5(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=1;j<dimension;j++){
            value[i][j] = sin(6*PI*value[i][0] + (j*PI)/dimension);
        }
    }

    fstream fpareto;
    string fparetoname = "UF5_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, hj1, sum2, yj2, hj2, N, E;
        
        sum1 = 0.0;
        count1 = 0;
        N = 10.0; E = 0.1;
        for(j = 2; j <= dimension; j++) {
            yj1 = value[i][j-1]-sin(6.0*PI*value[i][0]+j*PI/dimension);
            hj1 = 2.0*yj1*yj1 - cos(4.0*PI*yj1) + 1.0;
            if (j % 2 == 1) {
                sum1  += hj1;
                count1++;
            }
        }
        hj1 = (0.5/N + E)*fabs(sin(2.0*N*PI*value[i][0]));
        pareto[i][0] =  (value[i][0] + hj1 + 2.0*sum1 / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            yj2 = value[i][j-1]-sin(6.0*PI*value[i][0]+j*PI/dimension);
            hj2 = 2.0*yj2*yj2 - cos(4.0*PI*yj2) + 1.0;
            if (j % 2 == 0) {
                sum2  += hj2;
                count2++;
            } 
        }
        hj2 = (0.5/N + E)*fabs(sin(2.0*N*PI*value[i][0]));
        pareto[i][1] = (1.0 - value[i][0] + hj2 + 2.0*sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF6(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=1;j<dimension;j++){
            value[i][j] = sin(6*PI*value[i][0] + (j*PI)/dimension);
        }
    }

    fstream fpareto;
    string fparetoname = "UF6_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, sum2, yj2;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++){
            yj1 = value[i][j-1] - sin(6.0*PI*value[i][0] + j*PI/dimension);
            yj1 = yj1 * yj1;
            if(j % 2 == 1) {
                sum1 += yj1;
                count1++;
            } 
        }
        pareto[i][0] = (value[i][0] + 2.0 * sum1 / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++){
            yj2 = value[i][j-1] - sin(6.0*PI*value[i][0] + j*PI/dimension);
            yj2 = yj2 * yj2;
            if(j % 2 == 0) {
                sum2 += yj2;
                count2++;
            } 
        }
        
        pareto[i][1] = (1.0 - sqrt(value[i][0]) + 2.0 * sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}

void opt::UF7(){
    double min = 0;
    double max = 1;
    int dimension = 30;
    int objective = 2;
    
    //value assign
    value.assign(opt_num, vector<double>(dimension, min));

    for(int i=1;i<opt_num;i++){       
        value[i][0] = value[i-1][0] + (max-min)/(opt_num-1);
    }
    for(int i=0;i<opt_num;i++){
        for(int j=2;j<=dimension;j++){
            value[i][j-1] = sin(6.0*PI*value[i][0]+j*PI/dimension);
        }
    }

    fstream fpareto;
    string fparetoname = "UF7_opt.txt";
    fpareto.open(fparetoname, ios::out);
    if(!fpareto){
        cout<<"Fail to open file: "<<fparetoname<<endl;
    }
    
    //pareto front
    pareto.assign(opt_num, vector<double>(objective, 0));
    for(int i=0;i<opt_num;i++){
        unsigned int j, count1, count2;
        double sum1, yj1, sum2, yj2;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            yj1 = value[i][j-1] - sin(6.0*PI*value[i][0]+j*PI/dimension); 
            if (j % 2 == 1) {
                sum1  += yj1*yj1;
                count1++;
            }
        }
        yj1 = pow(value[i][0],0.2);
        pareto[i][0] = (yj1 + 2.0*sum1 / (double)count1);

        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            yj2 = value[i][j-1] - sin(6.0*PI*value[i][0]+j*PI/dimension);
            if (j % 2 == 0) {
                sum2  += yj2*yj2;
                count2++;
            } 
        }
        yj2 = pow(value[i][0],0.2);
        pareto[i][1] = (1.0 - yj2 + 2.0*sum2 / (double)count2);
        
        fpareto <<pareto[i][0]<<" "<<pareto[i][1]<<endl;

    }
}
