#!/bin/bash
REQWORKERS=100
WAITWORKERS=60
#RUNNUMBERS=" 246994"
RUNNUMBERS=$(cat runList.txt)
STATUSFILE=status_runs.txt
vafctl start
vafreq $REQWORKERS
rm -f $STATUSFILE
for RUN in $RUNNUMBERS; do
  vafwait $WAITWORKERS
  root -l -b -q "RunPoD.C($RUN)"
  RUNPAD=$(printf "%09d" $RUN)
  if [[ ! -e AnalysisResults_run$RUNPAD.root ]]; then
    echo "ATTENTION: Problem when running run $RUN. Cleaning up PROOF!"
    echo "$RUNPAD: BAD" >> $STATUSFILE
    vafctl stop
    vafctl start
    vafreq $REQWORKERS
  else
    echo "$RUNPAD: OK" >> $STATUSFILE
  fi
done
vafctl stop
