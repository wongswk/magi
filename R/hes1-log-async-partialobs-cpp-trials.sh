while true; do 
  Rscript R/hes1-log-async-partialobs-cpp-trials.R &> hes1-log-async-partialobs-cpp-trials-$(date +"%m-%d-%Y-%T").log
done
