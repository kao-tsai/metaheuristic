
#!bin/bash
clear
make clean
make

run=10
iter=100
cross_p=0.8
mut_p=0.2
pop_num=500
#./exe run iter dataset
#gnuplot -c exename filename dataname x1 x2 y1 y2 datafile

# ./search se $run 100 100 30 100 SCH

#  > pareto/FON/FON_se_0.txt
#250 100 8 2 6 4 good choice
#250 100 4 4 6 2 good choice
# ./search se $run 500 100 4 4 6 2 FON 0.9 0.033333
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_FON.png FON 0 1 0 1 pareto/FON/se/FON_se_0.txt opt/FON_opt.txt

# # # 250 500 0.85 0.4
# ./search se $run 500 100 4 4 6 2 SCH 0.9 0.033333
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_SCH.png SCH 0 4 0 4 pareto/SCH/se/SCH_se_0.txt opt/SCH_opt.txt

# ./search se $run 500 100 4 4 6 2 KUR 0.9 0.033333
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_KUR.png KUR -20 -14 -12 2 pareto/KUR/se/KUR_se_0.txt opt/KUR_opt.txt

# ./search se $run 500 100 4 4 6 2 POL 0.9 0.033333
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_POL.png POL 0 6 0 30 pareto/POL/se/POL_se_0.txt opt/POL_opt.txt

# ./search se $run 500 100 4 4 6 2 ZDT1 0.9 0.001
#  gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/se/ZDT1_se_0.txt opt/ZDT1_opt.txt

#  ./search se $run 500 100 4 4 6 2 ZDT2 0.9 0.001
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/se/ZDT2_se_0.txt opt/ZDT2_opt.txt

./search se $run 500 100 4 4 6 2 ZDT3 0.9 0.01
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/se/ZDT3_se_0.txt opt/ZDT3_opt.txt

./search se $run 500 100 4 4 6 2 ZDT4 0.9 0.01
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_ZDT4.png ZDT4 0 1 0 5 pareto/ZDT4/se/ZDT4_se_0.txt opt/ZDT4_opt.txt

./search se $run 500 100 4 4 6 2 ZDT6 0.9 0.0005
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/se/ZDT6_se_0.txt opt/ZDT6_opt.txt

./search se $run 500 100 4 4 6 2 UF1 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/se/UF1_se_0.txt opt/UF1_opt.txt

./search se $run 500 100 4 4 6 2 UF2 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF2.png UF2 0 1.2 0 1.2 pareto/UF2/se/UF2_se_0.txt opt/UF2_opt.txt

./search se $run 500 100 4 4 6 2 UF3 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF3.png UF3 0 1.2 0 1.2 pareto/UF3/se/UF3_se_0.txt opt/UF3_opt.txt

./search se $run 500 100 4 4 6 2 UF4 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF4.png UF4 0 1.2 0 1.2 pareto/UF4/se/UF4_se_0.txt opt/UF4_opt.txt

./search se $run 500 100 4 4 6 2 UF5 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF5.png UF5 0 1.2 0 1.2 pareto/UF5/se/UF5_se_0.txt opt/UF5_opt.txt

./search se $run 500 100 4 4 6 2 UF6 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF6.png UF6 0 1.2 0 1.2 pareto/UF6/se/UF6_se_0.txt opt/UF6_opt.txt

./search se $run 500 100 4 4 6 2 UF7 0.9 0.001
gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF7.png UF7 0 1.2 0 1.2 pareto/UF7/se/UF7_se_0.txt opt/UF7_opt.txt

# ./search se $run 500 100 4 4 6 2 UF8 0.9 0.001
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF8.png UF8 0 1.2 0 1.2 pareto/UF8/se/UF8_se_0.txt opt/UF8_opt.txt

# ./search se $run 500 100 4 4 6 2 UF9 0.9 0.001
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF9.png UF9 0 1.2 0 1.2 pareto/UF9/se/UF9_se_0.txt opt/UF9_opt.txt

# ./search se $run 500 100 4 4 6 2 UF10 0.9 0.001
# gnuplot -c gnuplot_se/test.gp gnuplot_se/SE_UF10.png UF10 0 1.2 0 1.2 pareto/UF10/se/UF10_se_0.txt opt/UF10_opt.txt


