#!/bin/bash
# worker: pops the first line of queue.txt (under a lock) and runs the corresponding rigorous prover
cd "$(dirname "$0")"
PY=/tmp/claude-0/-home-user-GapYas/76d3ef7a-1fe3-5a67-bc9e-3f6706c7d572/scratchpad/v/bin/python
while true; do
  job=$(flock queue.lock bash -c 'l=$(head -n1 queue.txt); [ -n "$l" ] && sed -i 1d queue.txt; echo "$l"')
  [ -z "$job" ] && exit 0
  set -- $job
  echo "$(date +%H:%M) start $job" >> queue.log
  if [ "$1" = rect ]; then $PY -u prove_rect_rig.py $2 $3 > r_rect$2.out 2>&1
  else $PY -u prove_min_rig.py $2 $1 > r_$1$2.out 2>&1; fi
  echo "$(date +%H:%M) end $job: $(grep -h '^DONE\|STALL\|Error' r_*$2.out | tail -n1 | cut -c1-80)" >> queue.log
done
