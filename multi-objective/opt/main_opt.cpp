#include "opt.h"

int main(int argc, char** argv) {

    int num = atoi(argv[1]);

    opt opt(num);
    opt.SCH();
    opt.FON();
    opt.KUR();
    opt.POL();
    opt.ZDT1();
    opt.ZDT2();
    opt.ZDT3();
    opt.ZDT4();
    opt.ZDT6();

    opt.UF1();
    opt.UF2();
    opt.UF3();
    opt.UF4();
    opt.UF7();

    /*
    opt.UF5();
    opt.UF6();
    */
}