#!/bin/bash

A=(1 2 4 8) # num motifs
B=(8 12 16 20) # motif length
C=(16 24 32 40) # num sequences
D=(256 512 1024 2048) # sequence length
LOG_FILE="serial_times.txt"
> $LOG_FILE  

for a in "${A[@]}"; do
    echo "Running for a=$a..." | tee -a $LOG_FILE
    total_time=0

    for i in {1..10}; do
        start_time=$(date +%s.%N)
        ./serial $a 16 16 $(512 - $(16 * $a))
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
        ./serial 2 $b 16 $(512 - $(2 * $a))
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
        ./serial 2 16 $c $(512 - 32)
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
        ./serial 2 16 16 $($d - 32)
        end_time=$(date +%s.%N)

        elapsed_time=$(echo "$end_time - $start_time" | bc)
        total_time=$(echo "$total_time + $elapsed_time" | bc)
    done

    average_time=$(echo "$total_time / 10" | bc -l)
    echo "Average time for c=$c: $average_time seconds" | tee -a $LOG_FILE
    echo | tee -a $LOG_FILE
done

echo "timing done"