make clean
make
# 3群4維
# ./search goaga 30 1024 "iris.txt" 20 3 0.8 0.2 1 0.00001>all_result/iris/goaga_iris_result.txt
# ./search hoo 30 500 30
#  > SE_based_HHO_result.txt
# ./search ga 30 1024 "iris.txt" 20 3 0.8 0.2>all_result/iris/ga_iris_result.txt

#./search 1.algo 2.runs 3.iteration 4.dimension 5.Region 6.searcher 7.sample 8.player
#較好 ./search se 30 500 30 2 2 20 4
# ./search se 30 200 30 1 8 10 2
# gnuplot -c gnuplot/test.gp gnuplot/compare.svg algo 0 257000 0 0.1 result/function/HHO.txt result/function/SEHHO.txt
# ./search tlbo 1 1000 30 30 Ackley

# MFO
# ./search mfo 30 10000 30 30 Ackley
# ./search mfo 30 10000 30 30 Rastrigin
./search mfo 30 10000 30 30 Sphere
# ./search mfo 30 10000 30 30 Rosenbrock
# ./search mfo 30 10000 30 30 Michalewicz
# ./search mfo 30 10000 30 30 Griewank10
# ./search mfo 30 10000 30 2 Schaffer2
# ./search mfo 30 10000 30 30 Schwefel
# ./search mfo 30 10000 30 2 Bohachevsky1
# ./search mfo 30 10000 30 30 Sum_Square
# ./search mfo 30 10000 30 2 Booth
# ./search mfo 30 10000 30 30 Zakharov
# ./search mfo 30 10000 30 2 Three_Hump_Camel
# ./search mfo 30 10000 30 2 De_Jong5
# ./search mfo 30 10000 30 2 Beale
# ./search mfo 30 10000 30 28 Powell

# MFODE
# ./search mfode 30 10000 30 30 Ackley
# ./search mfode 30 10000 30 30 Rastrigin
# ./search mfode 30 10000 30 30 Sphere
# ./search mfode 30 10000 30 30 Rosenbrock
# ./search mfode 30 10000 30 30 Michalewicz
# ./search mfode 30 10000 30 30 Griewank10
# ./search mfode 30 10000 30 2 Schaffer2
# ./search mfode 30 10000 30 30 Schwefel
# ./search mfode 30 10000 30 2 Bohachevsky1
# ./search mfode 30 10000 30 30 Sum_Square
# ./search mfode 30 10000 30 2 Booth
# ./search mfode 30 10000 30 30 Zakharov
# ./search mfode 30 10000 30 2 Three_Hump_Camel
# ./search mfode 30 10000 30 2 De_Jong5
# ./search mfode 30 10000 30 2 Beale
# ./search mfode 30 10000 30 28 Powell

#DE
# ./search de 30 10000 30 30 Ackley
# ./search de 30 10000 30 30 Rastrigin
# ./search de 30 10000 30 30 Sphere
# ./search de 30 10000 30 30 Rosenbrock
# ./search de 30 10000 30 30 Michalewicz
# ./search de 30 10000 30 30 Griewank10
# ./search de 30 10000 30 2 Schaffer2
# ./search de 30 10000 30 30 Schwefel
# ./search de 30 10000 30 2 Bohachevsky1
# ./search de 30 10000 30 30 Sum_Square
# ./search de 30 10000 30 2 Booth
# ./search de 30 10000 30 30 Zakharov
# ./search de 30 10000 30 2 Three_Hump_Camel
# ./search de 30 10000 30 2 De_Jong5
# ./search de 30 10000 30 2 Beale
# ./search de 30 10000 30 28 Powell

#PSO
# ./search pso 30 3000 30 30 Ackley
# ./search pso 30 3000 30 30 Rastrigin
# ./search pso 30 3000 30 30 Sphere
# ./search pso 30 3000 30 30 Rosenbrock
# ./search pso 30 3000 30 30 Michalewicz

#GA
# ./search ga 30 3000 30 30 0.9 0.01 Ackley
# ./search ga 30 3000 30 30 0.9 0.01 Rastrigin
# ./search ga 30 3000 30 30 0.9 0.01 Sphere
# ./search ga 30 3000 30 30 0.9 0.01 Rosenbrock
# ./search ga 30 3000 30 30 0.9 0.01 Michalewicz


# gnuplot -c gnuplot/test.gp gnuplot/Ackley_plot.svg Ackley 1 3000 0 20 result/mfo/Ackley.txt result/pso/Ackley.txt result/ga/Ackley.txt

# gnuplot -c gnuplot/test.gp gnuplot/Rastrigin_plot.svg Rastrigin 30 3000 40 400 result/mfo/Rastrigin.txt result/pso/Rastrigin.txt result/ga/Rastrigin.txt

# gnuplot -c gnuplot/test.gp gnuplot/Sphere_plot.svg Sphere 4 3000 0 40 result/mfo/Sphere.txt result/pso/Sphere.txt result/ga/Sphere.txt
# gnuplot -c gnuplot/test.gp gnuplot/Sphere_plot_b.svg Sphere 0.1 3000 0 1 result/mfo/Sphere.txt result/pso/Sphere.txt result/ga/Sphere.txt

# gnuplot -c gnuplot/test.gp gnuplot/Rosenbrock_plot.svg Rosenbrock 10000 3000 0 100000 result/mfo/Rosenbrock.txt result/pso/Rosenbrock.txt result/ga/Rosenbrock.txt
# gnuplot -c gnuplot/test.gp gnuplot/Michalewicz_plot.svg Michalewicz 1 3000 -29 -4 result/mfo/Michalewicz.txt result/pso/Michalewicz.txt result/ga/Michalewicz.txt

