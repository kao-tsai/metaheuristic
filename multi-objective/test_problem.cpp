#include "test_problem.h"
#define PI  3.1415926535897932384626433832795
using namespace std;

test_problem::test_problem(string func_name){
    if(func_name=="SCH"){
        idimension = 1;
        iobj_num=2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = -1000;
            ub[i] = 1000;
            xbits[i]=30;
        }
    }
    else if(func_name=="FON"){
        idimension = 3;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = -4;
            ub[i] = 4;
            xbits[i]=30;
        }
    }
    else if(func_name=="POL"){
        idimension = 2;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = -M_PI;
            ub[i] = M_PI;
            xbits[i]=30;
        }
    }
    else if(func_name=="KUR"){
        idimension = 3;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = -5;
            ub[i] = 5;
            xbits[i]=30;
        }
    }
    else if(func_name=="ZDT1"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = 0;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="ZDT2"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = 0;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="ZDT3"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = 0;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="ZDT4"){
        idimension = 10;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0]=0;
        ub[0]=1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -5;
            ub[i] = 5;
            xbits[i]=30;
        }
    }
    else if(func_name=="ZDT6"){
        idimension = 10;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = 0;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF1"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -1;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF2"){
        
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -1;
            ub[i] = 1;
            xbits[i]=30;
        }
        
    }
    else if(func_name=="UF3"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        for(int i=0;i<idimension;i++){
            lb[i] = 0;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF4"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -2;
            ub[i] = 2;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF5"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -1;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF6"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -1;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF7"){
        idimension = 30;
        iobj_num = 2;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = 0;
        ub[0] = 1;
        xbits[0]=30;
        for(int i=1;i<idimension;i++){
            lb[i] = -1;
            ub[i] = 1;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF8"){
        idimension = 30;
        iobj_num = 3;

        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = lb[1] = 0;
        ub[0] = ub[1] = 1;
        xbits[0]=xbits[1]=30;
        for(int i=2;i<idimension;i++){
            lb[i] = -2;
            ub[i] = 2;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF9"){
        idimension = 30;
        iobj_num = 3;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = lb[1] = 0;
        ub[0] = ub[1] = 1;
        xbits[0]=xbits[1]=30;
        for(int i=2;i<idimension;i++){
            lb[i] = -2;
            ub[i] = 2;
            xbits[i]=30;
        }
    }
    else if(func_name=="UF10"){
        idimension = 30;
        iobj_num = 3;
        ub.assign(idimension, 0);
        lb.assign(idimension, 0);
        xbits.assign(idimension,0);
        lb[0] = lb[1] = 0;
        ub[0] = ub[1] = 1;
        xbits[0]=xbits[1]=30;
        for(int i=2;i<idimension;i++){
            lb[i] = -2;
            ub[i] = 2;
            xbits[i]=30;
        }
    }
}
/*
double mofunc::mofunc_cal(vector<double> value, vector<double> &fit){
    if(func_name=="SCH"){
        double sum1 = 0;
        double sum2 = 0;
        
        sum1 = pow(value[0],2);
        fit[0] = sum1;
        
        sum2 = pow((value[0]-2),2);
        fit[1] = sum2;
    }
    else if(func_name=="FON"){
        double sum1 = 0;
        double sum2 = 0;
        
        for(int j=0;j<dimension;j++){
            sum1 += pow(value[j]-1/sqrt(3),2);
        }
        sum1 *= -1;
        sum1 = 1-exp(sum1);
        fit[0] = sum1;
            
        for(int j=0;j<dimension;j++){
            sum2 += pow(value[j]+1/sqrt(3),2);
        }
        sum2 *= -1;
        sum2 = 1-exp(sum2);
        fit[1] = sum2;
    }
    else if(func_name=="POL"){
        double sum1 = 0;
        double sum2 = 0;
        double A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
        double A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
        double B1 = 0.5*sin(value[0])-2*cos(value[0])+sin(value[1])-1.5*cos(value[1]);
        double B2 = 1.5*sin(value[0])-cos(value[0])+2*sin(value[1])-0.5*cos(value[1]);

        sum1 = 1+pow((A1-B1),2)+pow((A2-B2),2);
        fit[0] = sum1;
        
        sum2 = pow((value[0]+3),2)+pow((value[1]+1),2);
        fit[1] = sum2;
    }
    else if(func_name=="KUR"){
        double sum1 = 0;
        double sum2 = 0;
        
        for(int j=0;j<dimension-1;j++){
            sum1 += -10*exp(-0.2*sqrt(pow(value[j],2)+pow(value[j+1],2)));
        }
        fit[0] = sum1;
        
        for(int j=0;j<dimension;j++){
            sum2 += pow(abs(value[j]),0.8) + 5*sin(pow(value[j],3));
        }
        fit[1] = sum2;
    }
    else if(func_name=="ZDT1"){
        double sum1 = 0;
        double sum2 = 0;
        double g = 0;

        sum1 = value[0];
        fit[0] = sum1;

        for(int j=1;j<dimension;j++){
            g += value[j];
        }
        g = 1 + 9 * g / (dimension-1);
        
        sum2 = g*(1-sqrt(value[0]/g));
        fit[1] = sum2;
    }
    else if(func_name=="ZDT2"){
        double sum1 = 0;
        double sum2 = 0;
        double g = 0;

        sum1 = value[0];
        fit[0] = sum1;
        
        for(int j=1;j<dimension;j++){
            g += value[j];
        }
        g = 1 + (9 * g / (dimension-1));
        
        sum2 = g*(1-pow(value[0]/g, 2));
        fit[1] = sum2;
    }
    else if(func_name=="ZDT3"){
        double sum1 = 0;
        double sum2 = 0;
        double g = 0;

        sum1 = value[0];
        fit[0] = sum1;

        for(int j=1;j<dimension;j++){
            g += value[j];
        }
        g = 1 + 9 * g / (dimension-1);

        sum2 = g * (1 - sqrt(value[0]/g) - (value[0]/g)*sin(10*M_PI*value[0]));
        fit[1] = sum2;
    }
    else if(func_name=="ZDT4"){
        double sum1 = 0;
        double sum2 = 0;
        double g = 0;

        sum1 = value[0];
        fit[0] = sum1;

        for(int j=1;j<dimension;j++){
            g += pow(value[j],2) - 10*cos(4*M_PI*value[j]);
        }
        g = 1 + 10 * (dimension-1) + g;
            
        sum2 = g*(1-sqrt(value[0]/g));
        fit[1] = sum2;
    }
    else if(func_name=="ZDT6"){
        double sum1 = 0;
        double sum2 = 0;
        double g = 0;

        sum1 = 1 - exp(-4*value[0]) * pow(sin(6*M_PI*value[0]),6);
        fit[0] = sum1;
        
        for(int j=1;j<dimension;j++){
            g += value[j];
        }
        g = 1 + 9 * pow(g / (dimension-1),0.25);
            
        sum2 = g * (1 - pow(sum1/g,2));
        fit[1] = sum2;
    }
    else if(func_name=="UF1"){
        unsigned int j, count1, count2;
        double sum1, sum2, yj;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++){
            yj = value[j-1] - sin(6.0*PI*value[0] + j*PI/dimension);
            yj = yj * yj;
            if(j % 2 == 1) {
                sum1 += yj;
                count1++;
            } 
        }
        fit[0] = (value[0] + 2.0 * sum1 / (double)count1);
        
        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++){
            yj = value[j-1] - sin(6.0*PI*value[0] + j*PI/dimension);
            yj = yj * yj;
            if(j % 2 == 0) {
                sum2 += yj;
                count2++;
            } 
        }
        fit[1] = (1.0 - sqrt(value[0]) + 2.0 * sum2 / (double)count2);
    }
    else if(func_name=="UF2"){
        unsigned int j, count1, count2;
        double sum1, sum2, yj;
        
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            if(j % 2 == 1) {
                yj = value[j-1]-0.3*value[0]*(value[0]*cos(24.0*PI*value[0]+4.0*j*PI/dimension)+2.0)*cos(6.0*PI*value[0]+j*PI/dimension);
                sum1 += yj*yj;
                count1++;
            }
        }
        fit[0] = (value[0] + 2.0 * sum1 / (double)count1);
            
        sum2  = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            if(j % 2 == 0) {
                yj = value[j-1]-0.3*value[0]*(value[0]*cos(24.0*PI*value[0]+4.0*j*PI/dimension)+2.0)*sin(6.0*PI*value[0]+j*PI/dimension);
                sum2 += yj*yj;
                count2++;
            } 
        }
        fit[1] = (1.0 - sqrt(value[0]) + 2.0 * sum2 / (double)count2);
    }
    else if(func_name=="UF3"){
        unsigned int j, count1, count2;
        double sum1, prod1,sum2, prod2, yj, pj;
            
        sum1 = 0.0;
        count1 = 0;
        prod1 = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-pow(value[0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
            pj = cos(20.0*yj*PI/sqrt(j+0.0));
            if (j % 2 == 1) {
                sum1  += yj*yj;
                prod1 *= pj;
                count1++;
            }
        }
        fit[0] = (value[0] + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
            
        sum2 = 0.0;
        count2 = 0;
        prod2  = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-pow(value[0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
            pj = cos(20.0*yj*PI/sqrt(j+0.0));
            if (j % 2 == 0) {
                sum2  += yj*yj;
                prod2 *= pj;
                count2++;
            } 
            
        }
        fit[1] = (1.0 - sqrt(value[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
    }
    else if(func_name=="UF4"){
        unsigned int j, count1, count2;
        double sum1, sum2, yj, hj;
            
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
            if (j % 2 == 1) {
                sum1  += hj;
                count1++;
            }
        }
        fit[0] = (value[0] + 2.0*sum1 / (double)count1);
            
        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
            if (j % 2 == 0) {
                sum2  += hj;
                count2++;
            }
        }
        fit[1] = (1.0 - value[0]*value[0]	+ 2.0*sum2 / (double)count2);
    }
    else if(func_name=="UF5"){
        unsigned int j, count1, count2;
        double sum1, sum2, yj, hj, N, E;
            
        sum1 = 0.0;
        count1 = 0;
        N = 10.0; E = 0.1;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
            if (j % 2 == 1) {
                sum1  += hj;
                count1++;
            }
        }
        hj = (0.5/N + E)*fabs(sin(2.0*N*PI*value[0]));
        fit[0] = (value[0] + hj + 2.0*sum1 / (double)count1);
            
        sum2 = 0.0;
        count2 = 0;
        N = 10.0; E = 0.1;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
            if (j % 2 == 0) {
                sum2  += hj;
                count2++;
            } 
        }
        hj = (0.5/N + E)*fabs(sin(2.0*N*PI*value[0]));
        fit[1] = (1.0 - value[0] + hj + 2.0*sum2 / (double)count2);
    }
    else if(func_name=="UF6"){
        unsigned int j, count1, count2;
        double sum1, prod1, sum2, prod2, yj, hj, pj, N, E;
        N = 2.0; E = 0.1;

        sum1 = 0.0;
        count1 = 0;
        prod1 = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            pj = cos(20.0*yj*PI/sqrt(j+0.0));
            if (j % 2 == 1) {
                sum1  += yj*yj;
                prod1 *= pj;
                count1++;
            }
        }

        hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*value[0]);
        if(hj<0.0) hj = 0.0;
        fit[0] = (value[0] + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
        
        sum2 = 0.0;
        count2 = 0;
        prod2  = 1.0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1]-sin(6.0*PI*value[0]+j*PI/dimension);
            pj = cos(20.0*yj*PI/sqrt(j+0.0));
            if (j % 2 == 0) {
                sum2  += yj*yj;
                prod2 *= pj;
                count2++;
            } 
        }

        hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*value[0]);
        if(hj<0.0) hj = 0.0;
        fit[1] = (1.0 - value[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
    }
    else if(func_name=="UF7"){
        unsigned int j, count1, count2;
        double sum1, sum2, yj;
            
        sum1 = 0.0;
        count1 = 0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1] - sin(6.0*PI*value[0]+j*PI/dimension); 
            if (j % 2 == 1) {
                sum1  += yj*yj;
                count1++;
            }
        }
        yj = pow(value[0],0.2);
        fit[0] = (yj + 2.0*sum1 / (double)count1);
            
        sum2 = 0.0;
        count2 = 0;
        for(j = 2; j <= dimension; j++) {
            yj = value[j-1] - sin(6.0*PI*value[0]+j*PI/dimension);
            if (j % 2 == 0) {
                sum2  += yj*yj;
                count2++;
            } 
        }
        yj = pow(value[0],0.2);
        fit[1] = (1.0 - yj + 2.0*sum2 / (double)count2);
    }
    else if(func_name=="UF8"){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= dimension; j++) 
		{
			yj = value[j-1] - 2.0*value[1]*sin(2.0*PI*value[0]+j*PI/dimension);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		fit[0] = cos(0.5*PI*value[0])*cos(0.5*PI*value[1]) + 2.0*sum1 / (double)count1;
		fit[1] = cos(0.5*PI*value[0])*sin(0.5*PI*value[1]) + 2.0*sum2 / (double)count2;
		fit[2] = sin(0.5*PI*value[0]) + 2.0*sum3 / (double)count3;
    }
    else if(func_name=="UF9"){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, E;
		
		E = 0.1;
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= dimension; j++) 
		{
			yj = value[j-1] - 2.0*value[1]*sin(2.0*PI*value[0]+j*PI/dimension);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		yj = (1.0+E)*(1.0-4.0*(2.0*value[0]-1.0)*(2.0*value[0]-1.0));
		if(yj<0.0) yj = 0.0;
		fit[0] = 0.5*(yj + 2*value[0])*value[1] + 2.0*sum1 / (double)count1;
		fit[1] = 0.5*(yj - 2*value[0] + 2.0)*value[1] + 2.0*sum2 / (double)count2;
		fit[2] = 1.0 - value[1] + 2.0*sum3 / (double)count3;
    }
    else if(func_name=="UF10"){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= dimension; j++) 
		{
			yj = value[j-1] - 2.0*value[1]*sin(2.0*PI*value[0]+j*PI/dimension);
			hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
			if(j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += hj;
				count2++;
			}
			else
			{
				sum3  += hj;
				count3++;
			}
		}
		fit[0] = cos(0.5*PI*value[0])*cos(0.5*PI*value[1]) + 2.0*sum1 / (double)count1;
		fit[1] = cos(0.5*PI*value[0])*sin(0.5*PI*value[1]) + 2.0*sum2 / (double)count2;
		fit[2] = sin(0.5*PI*value[0]) + 2.0*sum3 / (double)count3;

    }

}*/
test_problem::population1 test_problem::SCH(population1 &all_sol){
    int obj_num=2;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    
    for(int i=0;i<all_sol.size();i++){
        all_sol_obj[i][0]=pow(all_sol[i][0],2.0);
        all_sol_obj[i][1]=pow(all_sol[i][0]-2.0,2.0);
    }
        
    return all_sol_obj;
}

test_problem::population1 test_problem::FON(population1 &all_sol){
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

test_problem::population1 test_problem::POL(population1 &all_sol){
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
test_problem::population1 test_problem::KUR(population1 &all_sol){
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
test_problem::population1 test_problem::ZDT1(population1 &all_sol){
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
test_problem::population1 test_problem::ZDT2(population1 &all_sol){
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
test_problem::population1 test_problem::ZDT3(population1 &all_sol){
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
test_problem::population1 test_problem::ZDT4(population1 &all_sol){
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
        // if(all_sol_obj[i][0]<0)
        //     cout<<g<<endl;
        f2 = g*(1-sqrt(all_sol[i][0]/g));
        all_sol_obj[i][1] = f2;
    }
    return all_sol_obj;
}

test_problem::population1 test_problem::ZDT6(population1 &all_sol){
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
test_problem::population1 test_problem::UF1(population1 &all_sol){
    int obj_num=2;
    int dimension=30;
    population1 all_sol_obj(all_sol.size(),solution1(obj_num));
    double yj,f1,f2,g;
    int count1,count2;
    for(int i=0;i<all_sol.size();i++){
        f1=0;
        f2=0;
        g=0;
        f1 = 0.0;
        count1 = 0;
        for(int j = 2; j <= dimension; j++){
            yj = all_sol[i][j-1] - sin(6.0*PI*all_sol[i][0] + j*PI/dimension);
            yj = yj * yj;
            if(j % 2 == 1) {
                f1 += yj;
                count1++;
            } 
        }
        all_sol_obj[i][0] = (all_sol[i][0] + 2.0 * f1 / (double)count1);
        
        f2 = 0.0;
        count2 = 0;
        for(int j = 2; j <= dimension; j++){
            yj = all_sol[i][j-1] - sin(6.0*PI*all_sol[i][0] + j*PI/dimension);
            yj = yj * yj;
            if(j % 2 == 0) {
                f2 += yj;
                count2++;
            } 
        }
        all_sol_obj[i][1] = (1.0 - sqrt(all_sol[i][0]) + 2.0 * f2 / (double)count2);
    }
    return all_sol_obj;
}

test_problem::population1 test_problem::UF2(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, sum2, yj;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            for(j = 2; j <= dimension; j++) {
                if(j % 2 == 1) {
                    yj = all_sol[i][j-1]-0.3*all_sol[i][0]*(all_sol[i][0]*cos(24.0*PI*all_sol[i][0]+4.0*j*PI/dimension)+2.0)*cos(6.0*PI*all_sol[i][0]+j*PI/dimension);
                    sum1 += yj*yj;
                    count1++;
                }
            }
            all_sol_obj[i][0] = (all_sol[i][0] + 2.0 * sum1 / (double)count1);
                
            sum2  = 0.0;
            count2 = 0;
            for(j = 2; j <= dimension; j++) {
                if(j % 2 == 0) {
                    yj = all_sol[i][j-1]-0.3*all_sol[i][0]*(all_sol[i][0]*cos(24.0*PI*all_sol[i][0]+4.0*j*PI/dimension)+2.0)*sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                    sum2 += yj*yj;
                    count2++;
                } 
            }
            all_sol_obj[i][1] = (1.0 - sqrt(all_sol[i][0]) + 2.0 * sum2 / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF3(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, prod1,sum2, prod2, yj, pj;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            prod1 = 1.0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-pow( all_sol[i][0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
                pj = cos(20.0*yj*PI/sqrt(j+0.0));
                if (j % 2 == 1) {
                    sum1  += yj*yj;
                    prod1 *= pj;
                    count1++;
                }
            }
            all_sol_obj[i][0] = (all_sol[i][0] + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
                
            sum2 = 0.0;
            count2 = 0;
            prod2  = 1.0;
            for(j = 2; j <= dimension; j++) {
                yj =  all_sol[i][j-1]-pow( all_sol[i][0],0.5*(1.0+3.0*(j-2.0)/(dimension-2.0)));
                pj = cos(20.0*yj*PI/sqrt(j+0.0));
                if (j % 2 == 0) {
                    sum2  += yj*yj;
                    prod2 *= pj;
                    count2++;
                } 
                
            }
            all_sol_obj[i][1] = (1.0 - sqrt(all_sol[i][0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF4(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, sum2, yj, hj;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI* all_sol[i][0]+(j*PI)/dimension);
                hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
                if (j % 2 == 1) {
                    sum1  += hj;
                    count1++;
                }
            }
            all_sol_obj[i][0] = (all_sol[i][0] + 2.0*sum1 / (double)count1);
                
            sum2 = 0.0;
            count2 = 0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
                if (j % 2 == 0) {
                    sum2  += hj;
                    count2++;
                }
            }
            all_sol_obj[j][1] = (1.0 -  all_sol[i][0]*all_sol[i][0]	+ 2.0*sum2 / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF5(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, sum2, yj, hj, N, E;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            N = 10.0; E = 0.1;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
                if (j % 2 == 1) {
                    sum1  += hj;
                    count1++;
                }
            }
            hj = (0.5/N + E)*fabs(sin(2.0*N*PI*all_sol[i][0]));
            all_sol_obj[i][0] = (all_sol[i][0] + hj + 2.0*sum1 / (double)count1);
                
            sum2 = 0.0;
            count2 = 0;
            N = 10.0; E = 0.1;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
                if (j % 2 == 0) {
                    sum2  += hj;
                    count2++;
                } 
            }
            hj = (0.5/N + E)*fabs(sin(2.0*N*PI*all_sol[i][0]));
            all_sol_obj[i][1] = (1.0 - all_sol[i][0] + hj + 2.0*sum2 / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF6(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, prod1, sum2, prod2, yj, hj, pj, N, E;
        N = 2.0; E = 0.1;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            prod1 = 1.0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                pj = cos(20.0*yj*PI/sqrt(j+0.0));
                if (j % 2 == 1) {
                    sum1  += yj*yj;
                    prod1 *= pj;
                    count1++;
                }
            }

            hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*all_sol[i][0]);
            if(hj<0.0) hj = 0.0;
            all_sol_obj[i][0] = (all_sol[i][0] + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1);
            
            sum2 = 0.0;
            count2 = 0;
            prod2  = 1.0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1]-sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                pj = cos(20.0*yj*PI/sqrt(j+0.0));
                if (j % 2 == 0) {
                    sum2  += yj*yj;
                    prod2 *= pj;
                    count2++;
                } 
            }

            hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*all_sol[i][0]);
            if(hj<0.0) hj = 0.0;
            all_sol_obj[i][1] = (1.0 - all_sol[i][0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF7(population1 &all_sol){
        unsigned int j, count1, count2;
        double sum1, sum2, yj;
        int obj_num=2;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1 = 0.0;
            count1 = 0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1] - sin(6.0*PI*all_sol[i][0]+j*PI/dimension); 
                if (j % 2 == 1) {
                    sum1  += yj*yj;
                    count1++;
                }
            }
            yj = pow(all_sol[i][0],0.2);
            all_sol_obj[i][0] = (yj + 2.0*sum1 / (double)count1);
                
            sum2 = 0.0;
            count2 = 0;
            for(j = 2; j <= dimension; j++) {
                yj = all_sol[i][j-1] - sin(6.0*PI*all_sol[i][0]+j*PI/dimension);
                if (j % 2 == 0) {
                    sum2  += yj*yj;
                    count2++;
                } 
            }
            yj = pow(all_sol[i][0],0.2);
            all_sol_obj[i][1] = (1.0 - yj + 2.0*sum2 / (double)count2);
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF8(population1 &all_sol){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;
		int obj_num=3;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1   = sum2   = sum3   = 0.0;
            count1 = count2 = count3 = 0;
            for(j = 3; j <= dimension; j++) 
            {
                yj = all_sol[i][j-1] - 2.0*all_sol[i][1]*sin(2.0*PI*all_sol[i][0]+j*PI/dimension);
                if(j % 3 == 1) 
                {
                    sum1  += yj*yj;
                    count1++;
                } 
                else if(j % 3 == 2) 
                {
                    sum2  += yj*yj;
                    count2++;
                }
                else
                {
                    sum3  += yj*yj;
                    count3++;
                }
            }
            all_sol_obj[i][0] = cos(0.5*PI*all_sol[i][0])*cos(0.5*PI*all_sol[i][1]) + 2.0*sum1 / (double)count1;
            all_sol_obj[i][1] = cos(0.5*PI*all_sol[i][0])*sin(0.5*PI*all_sol[i][1]) + 2.0*sum2 / (double)count2;
            all_sol_obj[i][2] = sin(0.5*PI*all_sol[i][0]) + 2.0*sum3 / (double)count3;
         }
         return all_sol_obj;
    }
    test_problem::population1 test_problem::UF9(population1 &all_sol){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, E;
		E = 0.1;
        int obj_num=3;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));

        for(int i=0;i<all_sol.size();i++){
            sum1   = sum2   = sum3   = 0.0;
            count1 = count2 = count3 = 0;
            for(j = 3; j <= dimension; j++) 
            {
                yj = all_sol[i][j-1] - 2.0* all_sol[i][1]*sin(2.0*PI* all_sol[i][0]+j*PI/dimension);
                if(j % 3 == 1) 
                {
                    sum1  += yj*yj;
                    count1++;
                } 
                else if(j % 3 == 2) 
                {
                    sum2  += yj*yj;
                    count2++;
                }
                else
                {
                    sum3  += yj*yj;
                    count3++;
                }
            }
            yj = (1.0+E)*(1.0-4.0*(2.0*all_sol[i][0]-1.0)*(2.0*all_sol[i][0]-1.0));
            if(yj<0.0) yj = 0.0;
            all_sol_obj[i][0] = 0.5*(yj + 2* all_sol[i][0])* all_sol[i][1] + 2.0*sum1 / (double)count1;
            all_sol_obj[i][1] = 0.5*(yj - 2* all_sol[i][0] + 2.0)* all_sol[i][1] + 2.0*sum2 / (double)count2;
            all_sol_obj[i][2] = 1.0 - all_sol[i][1] + 2.0*sum3 / (double)count3;
        }
        return all_sol_obj;
    }
    test_problem::population1 test_problem::UF10(population1 &all_sol){
        unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj;
		int obj_num=3;
        int dimension=30;
        population1 all_sol_obj(all_sol.size(),solution1(obj_num));
        for(int i=0;i<all_sol.size();i++){
            sum1   = sum2   = sum3   = 0.0;
            count1 = count2 = count3 = 0;
            for(j = 3; j <= dimension; j++) 
            {
                yj = all_sol[i][j-1] - 2.0*all_sol[i][1]*sin(2.0*PI*all_sol[i][0]+j*PI/dimension);
                hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
                if(j % 3 == 1) 
                {
                    sum1  += hj;
                    count1++;
                } 
                else if(j % 3 == 2) 
                {
                    sum2  += hj;
                    count2++;
                }
                else
                {
                    sum3  += hj;
                    count3++;
                }
            }
            all_sol_obj[i][0] = cos(0.5*PI*all_sol[i][0])*cos(0.5*PI*all_sol[i][1]) + 2.0*sum1 / (double)count1;
            all_sol_obj[i][1] = cos(0.5*PI*all_sol[i][0])*sin(0.5*PI*all_sol[i][1]) + 2.0*sum2 / (double)count2;
            all_sol_obj[i][2] = sin(0.5*PI*all_sol[i][0]) + 2.0*sum3 / (double)count3;
        }
        return all_sol_obj;
    }