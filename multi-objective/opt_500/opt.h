#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <sstream>

using namespace std;

#define PI  3.1415926535897932384626433832795

class opt{
    private:
        int opt_num;
        vector<vector<double> > value;
        vector<vector<double> > pareto;
    public:
        opt(int num);
        void value_assign();
        void SCH();
        void FON();
        void KUR();
        void POL();
        void ZDT1();
        void ZDT2();
        void ZDT3();
        void ZDT4();
        void ZDT6();

        void UF1();
        void UF2();
        void UF3();
        void UF4();
        void UF5();
        void UF6();
        void UF7();

};