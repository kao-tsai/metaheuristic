#ifndef __test_problem_H_INCLUDED__
#define __test_problem_H_INCLUDED__
#include <fstream>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <climits>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define PI  3.1415926535897932384626433832795
// #define PI  3.14159
class test_problem{
    public:
        test_problem(){};
        test_problem(string,int);
        typedef vector<double> solution1;
        typedef vector<solution1> population1;
        
        solution1 ub;
        solution1 lb;
        vector<int> xbits;
        double Ackley(solution1&);
        double Rastrigin(solution1&);
        double Sphere(solution1&);
        double Rosenbrock(solution1&);
        double Michalewicz(solution1&);
        double F15(solution1&);
        double F16(solution1&);
        int idim;

}; 
#endif