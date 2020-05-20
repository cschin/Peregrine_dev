#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate peregrine
if [ $1 == "test" ]; then 
  cd /opt/test
  bash /opt/test/run_test.sh
else  
  python /opt/generate_cns.py $@
fi
