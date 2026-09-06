#!/bin/bash
cd /Users/jrines/stanford_gp/research/lakes/2026/persistent_ids
while ! grep -q "batch done" out/08_batch_29_45.log; do sleep 300; done
n=$(ls out/08_s2counts_29_45_*_meta.json 2>/dev/null | wc -l); echo "$(TZ=America/Los_Angeles date '+%H:%M') batch done, $n seasons downloaded; running 08e" >> out/08_batch_29_45.log
TILE=29_45 TZ=America/Los_Angeles nice -n 15 $(cat .python_env) scripts/08e_sites_registry.py >> out/08_batch_29_45.log 2>&1
echo "$(TZ=America/Los_Angeles date '+%H:%M') 08e done" >> out/08_batch_29_45.log
