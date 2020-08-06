while true; do 
  Rscript R/hes1-log-async-partialobs-cpp-trials.R &> hes1-log-async-partialobs-cpp-trials-$(date +"%m-%d-%Y-%T").log
done

while true; do 
  Rscript R/hes1-log-fullobs-cpp.R

sleep 2400

for dummy in {1..12}; do
  for i in {1..60}; do
    Rscript R/hes1-log-async-partialobs-cpp-temper.R &
  done
  Rscript R/hes1-log-async-partialobs-cpp-temper.R
  sleep 300
done

for dummy in {1..12}; do
  for i in {1..60}; do
    Rscript  R/fn-model-cpp.R &
  done
  Rscript  R/fn-model-cpp.R
  sleep 60
done

iseed=1
for dummy in {1..12}; do
  for i in {1..60}; do
    Rscript R/ramsey-rerun-hes1-partialobs.R $iseed &
    ((iseed=iseed+1))
  done
  Rscript R/ramsey-rerun-hes1-partialobs.R $iseed
  ((iseed=iseed+1))
  sleep 300
done


sleep 10800
for i in {1..60}
do
  Rscript R/hes1-log-async-partialobs-cpp-trials.R &
done

for i in {193..256}
do
  echo "Number: $i"
  Rscript R/ramsey-rerun-hes1-partialobs.R $i
done

for dummy in {1..3}; do
  for i in {1..60}; do
    Rscript R/fn-temper.R &
  done
  Rscript R/fn-temper.R
  sleep 180
done

for dummy in {1..3}; do
  for i in {1..60}; do
    Rscript R/ptrans-temper.R &
  done
  Rscript R/ptrans-temper.R
  sleep 180
done
