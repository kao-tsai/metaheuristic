#ifndef __test_problem_H_INCLUDED__
#define __test_problem_H_INCLUDED__
#include <fstream>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <climits>
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

#define PI  3.1415926535897932384626433832795

class test_problem{
    public:
        test_problem(){};
        test_problem(string);
        typedef vector<double> solution1;
        typedef vector<solution1> population1;
        int idimension;
        int iobj_num;
        solution1 ub;
        solution1 lb;
        vector<int> xbits;
    
        population1 SCH(population1&);
        population1 FON(population1&);
        population1 KUR(population1&);
        population1 POL(population1&);
        population1 ZDT1(population1&);
        population1 ZDT2(population1&);
        population1 ZDT3(population1&);
        population1 ZDT4(population1&);
        population1 ZDT6(population1&);

        population1 UF1(population1&);
        population1 UF2(population1&);
        population1 UF3(population1&);
        population1 UF4(population1&);
        population1 UF5(population1&);
        population1 UF6(population1&);
        population1 UF7(population1&);
        population1 UF8(population1&);
        population1 UF9(population1&);
        population1 UF10(population1&);

}; 
#endif