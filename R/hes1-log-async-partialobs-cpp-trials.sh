while true; do 
  Rscript R/hes1-log-async-partialobs-cpp-trials.R &> hes1-log-async-partialobs-cpp-trials-$(date +"%m-%d-%Y-%T").log
done

while true; do 
  Rscript R/hes1-log-fullobs-cpp.R
done

for i in {193..256}
do
  echo "Number: $i"
  Rscript R/ramsey-rerun-hes1-partialobs.R $i
done
