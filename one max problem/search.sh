clear
make clean
make
#alpha fine tune
./search aco 30 1024 "ACO/testData.txt" 20 0.1 2 0.1 0
./search aco 30 1024 "ACO/testData.txt" 20 1 2 0.1 0
./search aco 30 1024 "ACO/testData.txt" 20 2 2 0.1 0
./search aco 30 1024 "ACO/testData.txt" 20 4 2 0.1 0
#beta fine tune
./search aco 30 1024 "ACO/testData.txt" 20 1 0.1 0.1 1
./search aco 30 1024 "ACO/testData.txt" 20 1 1 0.1 1
./search aco 30 1024 "ACO/testData.txt" 20 1 2 0.1 1
./search aco 30 1024 "ACO/testData.txt" 20 1 4 0.1 1
#rho fine tune
./search aco 30 1024 "ACO/testData.txt" 20 1 0.1 0.02 2
./search aco 30 1024 "ACO/testData.txt" 20 1 1 0.1 2
./search aco 30 1024 "ACO/testData.txt" 20 1 2 0.5 2
./search aco 30 1024 "ACO/testData.txt" 20 1 4 0.9 2

bash ACO/plot_each_para.sh
