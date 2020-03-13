for ((i=0; i<1000; i++));do
  jname=$((${i} + 1))
  for ((j=0; j<10; j++));do
     cname=$((${j} + 1))
     cd "/mnt/home/thrust2/zf1005/Matlab/GAFitCyclic/RunningFolder/$((${i} + 1))/$((${j} + 1))/" && chmod 777 a.out
  done
done
