while true; do 
  Rscript R/hes1-log-async-partialobs-cpp-trials.R &> hes1-log-async-partialobs-cpp-trials-$(date +"%m-%d-%Y-%T").log
done

while true; do 
  Rscript R/hes1-log-fullobs-cpp.R

sleep 2400
for i in {1..60}
do
  Rscript R/hes1-log-async-partialobs-cpp-temper.R &
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
