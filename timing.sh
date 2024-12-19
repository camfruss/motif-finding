#!/bin/bash

A=(1 4 16 64) # num motifs
B=(8 16 32 64) # motif length
C=(16 32 64 128) # num sequences
D=(128 512 2048 8192) # sequence length
LOG_FILE="serial_times.txt"
> $LOG_FILE  

for a in "${A[@]}"; do
    echo "Running for a=$a..." | tee -a $LOG_FILE
    total_time=0

    for i in {1..10}; do
        start_time=$(date +%s.%N)
        ./serial $a 16 64 512
        end_time=$(date +%s.%N)

        elapsed_time=$(echo "$end_time - $start_time" | bc)
        total_time=$(echo "$total_time + $elapsed_time" | bc)
    done

    average_time=$(echo "$total_time / 10" | bc -l)
    echo "Average time for a=$a: $average_time seconds" | tee -a $LOG_FILE
    echo | tee -a $LOG_FILE
done

for b in "${B[@]}"; do
    echo "Running for b=$b..." | tee -a $LOG_FILE
    total_time=0

    for i in {1..10}; do
        start_time=$(date +%s.%N)
        ./serial 4 $b 64 512
        end_time=$(date +%s.%N)

        elapsed_time=$(echo "$end_time - $start_time" | bc)
        total_time=$(echo "$total_time + $elapsed_time" | bc)
    done

    average_time=$(echo "$total_time / 10" | bc -l)
    echo "Average time for b=$b: $average_time seconds" | tee -a $LOG_FILE
    echo | tee -a $LOG_FILE
done

for c in "${C[@]}"; do
    echo "Running for c=$c..." | tee -a $LOG_FILE
    total_time=0

    for i in {1..10}; do
        start_time=$(date +%s.%N)
        ./serial 4 16 $c 512
        end_time=$(date +%s.%N)

        elapsed_time=$(echo "$end_time - $start_time" | bc)
        total_time=$(echo "$total_time + $elapsed_time" | bc)
    done

    average_time=$(echo "$total_time / 10" | bc -l)
    echo "Average time for c=$c: $average_time seconds" | tee -a $LOG_FILE
    echo | tee -a $LOG_FILE
done

for d in "${D[@]}"; do
    echo "Running for d=$d..." | tee -a $LOG_FILE
    total_time=0

    for i in {1..10}; do
        start_time=$(date +%s.%N)
        ./serial 4 16 64 $d
        end_time=$(date +%s.%N)

        elapsed_time=$(echo "$end_time - $start_time" | bc)
        total_time=$(echo "$total_time + $elapsed_time" | bc)
    done

    average_time=$(echo "$total_time / 10" | bc -l)
    echo "Average time for c=$c: $average_time seconds" | tee -a $LOG_FILE
    echo | tee -a $LOG_FILE
done

echo "timing done"