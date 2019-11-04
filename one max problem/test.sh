clear
make clean
make
./search aco 30 1024 "ACO/testData.txt" 20 1 2 0.9 4 > result_tsp.txt

bash ACO/plot_best.sh
