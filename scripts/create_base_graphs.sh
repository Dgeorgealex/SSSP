#!/bin/bash

for sizes in "2 6"; do
    size=( $sizes )
    coeff=${size[0]}
    pow=${size[1]}

    for alg in "gor"; do
        for ((iter=1; iter<=3; iter++)); do

            if [ ! -f "../data/graphs/big_aug_${alg}_${coeff}e${pow}_v${iter}.txt" ]; then
                bash big_graph_creator.sh $coeff $pow ${alg} 11

                # Delete temporary files
                rm ../data/graphs/big_aug_${alg}_${coeff}e${pow}_bare.txt
                rm ../data/graphs/big_aug_${alg}_${coeff}e${pow}_sorted.txt
                mv ../data/graphs/big_aug_${alg}_${coeff}e${pow}.txt ../data/graphs/big_aug_${alg}_${coeff}e${pow}_v${iter}.txt
            else
                echo "../data/graphs/big_aug_${alg}_${coeff}e${pow}_v${iter}.txt already exists, skipping..."
            fi
        done

    done
done