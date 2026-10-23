#!/bin/bash
# Run solver2 with many weight configs, saving per-pattern solutions to sols_large
cd /home/user/CCC3/level6
B=${1:-100}
while read -r w; do
  [ -z "$w" ] && continue
  echo "== $w" >> portfolio.log
  MAXD=250 SOLDIR=sols_large ./solver2 in_level-6_3-large.txt /dev/null $B 4 $w 2>>portfolio.log
  python3 -I merge.py in_level-6_3-large.txt sols_large merged.txt | tail -2 >> portfolio.log
done <<'CFG'
30 3 6 150 0.5 3
30 1 2 100 0.5 2
30 3 6 60 0.5 1
30 10 20 100 0.5 2
30 6 12 200 0.5 3
30 2 6 300 0.5 4
30 4 8 120 0.5 2
30 1 4 200 0.5 3
CFG
echo ALLDONE >> portfolio.log
