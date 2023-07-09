#set -x

test=./Executable/stest_float
test2=./Executable/stest2x2_float
test4=./Executable/stest4x4_float
test8=./Executable/stest8x8_float

echo $test
for N in 2000 4000 6000
do

    for D in 16 32 64
    do

	for p in 4 6 8 10 12 14
	do
	    echo $N $D $p | $test | grep GFLOPS | grep -v BLAS
    
	done
    done
done


echo $test2
for N in 1000 2000 3000 
do

    for D in 8 16 32
    do

	for p in 4 6 8 10 12 14
	do
	    echo $N $D $p | $test2 | grep GFLOPS | grep -v BLAS
    
	done
    done
done

#exit 
echo $test4
for N in  500 1000 2000
do

    for D in 1 8 16 
    do

	for p in 4 6 8 10 12 14
	do
	    echo $N $D $p | $test4 | grep GFLOPS | grep -v BLAS
    
	done
    done
done

echo $test8
for N in 250 500 1000 
do

    for D in 8 16
    do

	for p in 2 4 6 8 10 12 14 16
	do

	    echo $N $D $p | $test8 | grep GFLOPS | grep -v BLAS
    
	done
    done
done
