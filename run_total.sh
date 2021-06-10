#! /bin/bash

set -o xtrace

MEASUREMENTS=10
ITERATIONS=10
INITIAL_SIZE=16
INITIAL_NUM_THREADS=1

SIZE=$INITIAL_SIZE
THREADS=$INITIAL_NUM_THREADS

NAMES=('mandelbrot_seq')
NAMESPARALLEL=('mandelbrot_pth' 'mandelbrot_omp')

make
mkdir results_minimal

for NAME in ${NAMES[@]}; do
    mkdir results_minimal/$NAME

    for ((i=1; i<=$ITERATIONS; i++)); do
            for ((j=1; j<=$MEASUREMENTS; j++)); do
                ./$NAME -2.5 1.5 -2.0 2.0 $SIZE >> full_minimal.log 2>&1
                ./$NAME -0.8 -0.7 0.05 0.15 $SIZE >> seahorse_minimal.log 2>&1
                ./$NAME 0.175 0.375 -0.1 0.1 $SIZE >> elephant_minimal.log 2>&1
                ./$NAME -0.188 -0.012 0.554 0.754 $SIZE >> triple_spiral_minimal.log 2>&1
            done
            echo "\n" >> full_minimal.log
            echo "\n" >> seahorse_minimal.log
            echo "\n" >> elephant_minimal.log
            echo "\n" >> triple_spiral_minimal.log
            SIZE=$(($SIZE * 2))
    done

    SIZE=$INITIAL_SIZE

    mv *.log results_minimal/$NAME
    rm output.ppm
done


for NAME in ${NAMESPARALLEL[@]}; do
    mkdir results_minimal/$NAME
    for ((p=1; p<=6; p++)); do
        SIZE=$INITIAL_SIZE
        for ((i=1; i<=$ITERATIONS; i++)); do
                for ((j=1; j<=$MEASUREMENTS; j++)); do
                    ./$NAME -2.5 1.5 -2.0 2.0 $SIZE $THREADS >> full_minimal.log 2>&1
                    ./$NAME -0.8 -0.7 0.05 0.15 $SIZE $THREADS >> seahorse_minimal.log 2>&1
                    ./$NAME 0.175 0.375 -0.1 0.1 $SIZE $THREADS >> elephant_minimal.log 2>&1
                    ./$NAME -0.188 -0.012 0.554 0.754 $SIZE $THREADS >> triple_spiral_minimal.log 2>&1
                done
                echo "\n" >> full_minimal.log
                echo "\n" >> seahorse_minimal.log
                echo "\n" >> elephant_minimal.log
                echo "\n" >> triple_spiral_minimal.log
                SIZE=$(($SIZE * 2))
        done
        THREADS=$(($THREADS * 2))
    done

    SIZE=$INITIAL_SIZE
    THREADS=$INITIAL_NUM_THREADS

    mv *.log results_minimal/$NAME
    rm output.ppm
done