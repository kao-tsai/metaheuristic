make clean
make
# 3群4維
# ./search goaga 30 1024 "iris.txt" 20 3 0.8 0.2 1 0.00001>all_result/iris/goaga_iris_result.txt
./search goa2 1 1024 "iris.txt" 20 3 1 0.00001>all_result/iris/goa_iris_result.txt
# ./search ga 30 1024 "iris.txt" 20 3 0.8 0.2>all_result/iris/ga_iris_result.txt
#6群12維
# ./search goaga 1 1000 "vowel.txt" 20 6 0.8 0.2 1 0.00001>all_result/vowel/goaga_vowel_result.txt
# ./search goa 1 1000 "vowel.txt" 20 6 1 0.00001>all_result/vowel/goa_vowel_result.txt
# ./search ga 1 1000 "vowel.txt" 20 6 0.8 0.2>all_result/vowel/ga_vowel_result.txt
#28群8維
# ./search goaga 1 10000 "abalone.txt" 20 28 0.8 0.2 1 0.00001>all_result/abalone/goaga_abalone_result.txt
# ./search goa 1 10000 "abalone.txt" 20 28 1 0.00001>all_result/abalone/goa_abalone_result.txt
# ./search ga 1 10000 "abalone.txt" 20 28 0.8 0.2>all_result/abalone/ga_abalone_result.txt
#3群13維
# ./search goaga 30 1024 "wine.txt" 20 3 0.8 0.2 1 0.00001>all_result/wine/goaga_wine_result.txt
# ./search goa 30 1024 "wine.txt" 20 3 1 0.00001>all_result/wine/goa_wine_result.txt
# ./search ga 30 1024 "wine.txt" 20 3 0.8 0.2>all_result/wine/ga_wine_result.txt