#!/bin/bash
echo "Unzipping USA-road-d.W.gr.gz..."
gunzip ../data/graphs/USA-road-d.W.gr.gz data/graphs/USA-road-d.W.gr

if [ ! -f "../data/graphs/USA-road-d.W.txt" ]; then
    echo "Converting from dimacs..."
    python3 from_dimacs.py ../data/graphs/USA-road-d.W.gr ../data/graphs/USA-road-d.W.txt
else
    echo "../data/graphs/USA-road-d.W.txt already exists, skipping..."
fi

echo "Generating feasible potential..."
# Solve the positive instance
../build/Main "SSSP time GOR ../data/graphs/USA-road-d.W.txt 0 1"

#
#if [ ! -f "../data/graphs/USA-road-d.W_1.txt" ]; then
#    echo "../data/graphs/USA-road-d.W_1.txt does not exist, creating..."
#    python3 create_graph.py restr_from_pot ../data/graphs/USA-road-d.W.txt ../data/graphs/USA-road-d.W_result0.txt 1 1 > ../data/graphs/USA-road-d.W_1.txt
#else
#    echo "../data/graphs/USA-road-d.W_1.txt already exists, skipping..."
#fi

if [ ! -f "../data/graphs/USA-road-d.W_10.txt" ]; then
    echo "../data/graphs/USA-road-d.W_10.txt does not exist, creating..."
    python3 create_graph.py restr_from_pot ../data/graphs/USA-road-d.W.txt ../data/graphs/USA-road-d.W_result0.txt 10 1 > ../data/graphs/USA-road-d.W_10.txt
else
    echo "../data/graphs/USA-road-d.W_10.txt already exists, skipping..."
fi

if [ ! -f "../data/graphs/USA-road-d.W_100.txt" ]; then
    echo "../data/graphs/USA-road-d.W_100.txt does not exist, creating..."
    python3 create_graph.py restr_from_pot ../data/graphs/USA-road-d.W.txt ../data/graphs/USA-road-d.W_result0.txt 100 1 > ../data/graphs/USA-road-d.W_100.txt
else
    echo "../data/graphs/USA-road-d.W_100.txt already exists, skipping..."
fi

if [ ! -f "../data/graphs/USA-road-d.W_500.txt" ]; then
    echo "../data/graphs/USA-road-d.W_500.txt does not exist, creating..."
    python3 create_graph.py restr_from_pot ../data/graphs/USA-road-d.W.txt ../data/graphs/USA-road-d.W_result0.txt 500 1 > ../data/graphs/USA-road-d.W_500.txt
else
    echo "../data/graphs/USA-road-d.W_500.txt already exists, skipping..."
fi
