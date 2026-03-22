#!/bin/bash
set -euo pipefail

for f in ../data/graphs/*.base; do
  [ -e "$f" ] || continue
  n=$(basename "$f" .base)
  ../build/CreateGraph load_graph "$f" -t -c 1 > "../data/graphs/${n}_t_c1.txt"
done