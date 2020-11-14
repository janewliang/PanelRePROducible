#!/bin/bash                                                                               

# Define the stats types to be run                                                        
type="diagnostics_boot"

for a in $type ; do
   for i in {1..100} ; do
		index="$index $i"
   done
   elements=`echo $index | tr " " ","`
   #if the array element is non-zero submit the job                                       
   if [ -n "$elements" ] ; then
      echo "Submitting command: "
      echo "sbatch --array=$elements rscript_${a}.job"
      sbatch --array=$elements rscript_${a}.job
   fi
   #null the arrays for next type                                                         
   index=""
   elements=""
done; 
wait
