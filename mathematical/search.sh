make clean
make
# 3群4維
# ./search goaga 30 1024 "iris.txt" 20 3 0.8 0.2 1 0.00001>all_result/iris/goaga_iris_result.txt
# ./search hoo 30 3500 30
# ./search ga 30 1024 "iris.txt" 20 3 0.8 0.2>all_result/iris/ga_iris_result.txt

#./search 1.algo 2.runs 3.iteration 4.dimension 5.Region 6.searcher 7.sample 8.player
#較好 ./search se 30 500 30 2 2 20 4
# ./search se 30 1000 30 2 10 5 2
gnuplot -c gnuplot/test.gp gnuplot/compare.svg algo 0 257000 0 0.1 result/function/HHO.txt result/function/SEHHO.txt

#  > SE_based_HHO_result.txt
