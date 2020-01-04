clear
make clean
make
./search ga 30 1024 "GA/testData.txt" 100 0.8 0.2 ts pmx > PMX_tsp_result.txt
./search ga 30 1024 "GA/testData.txt" 100 0.8 0.2 ts cx > CX_tsp_result.txt
./search ga 30 1024 "GA/testData.txt" 100 0.8 0.2 ts ox > OX_tsp_result.txt

