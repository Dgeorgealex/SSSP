#!/bin/bash
set -euo pipefail

for f in ../data/graphs/*.base; do
  [ -e "$f" ] || continue
  n=$(basename "$f" .base)
  ../build/CreateGraph load_graph "$f" -c 2 -t > "../data/graphs/${n}_c2_0.txt"
done