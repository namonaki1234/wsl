#!/bin/sh
#PBS -l nodes=yato215:ppn=16
#PBS -N ABE_ALL
#PBS -o /dev/null
#PBS -j oe

MULTISHOT_MAX=10

OMP_NUM_THREADS=16; export OMP_NUM_THREADS
ulimit -s unlimited
cd $PBS_O_WORKDIR

#rm ./log/GridFlow/*.txt
#rm ./log/MPSini/*.txt
#rm ./log/MPS/*.txt
#rm ./log/ICM/*.txt

RESULTDIR="202401270050"
SEQNUM=7
#mkdir -p result/$RESULTDIR

(echo ==========================	) >  ./log.txt
(echo Computation Start			) >> ./log.txt
(echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt
(echo ==========================	) >> ./log.txt

(echo 8; echo 2) |
for i in `seq -w ${SEQNUM} ${MULTISHOT_MAX}`
do
(echo $i) > ./data/STEP.txt
(echo ${i}/${MULTISHOT_MAX}+++++++++++++++++++++++	) >> ./log.txt

 #GRID_FLOW
 if [ $i -eq 1 ]; then
  (echo Start Grid, Overset and InitialFlow	) >> ./log.txt
  (time ./temp/GridFlow/Grid.exe      ) &> ./log/GridFlow/log_step${i}_Grid.txt
  (time ./temp/GridFlow/Overset.exe   ) &> ./log/GridFlow/log_step${i}_Overset.txt
  (time ./temp/GridFlow/IniFlow.exe    ) &> ./log/GridFlow/log_step${i}_IniFlow.txt
  rsync -r ./result/GridFlow/grid ./result/${RESULTDIR}
  rsync -r ./result/GridFlow/overset ./result/${RESULTDIR}
  (echo Complete Grid, Overset and InitialFlow	) >> ./log.txt
  (echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt
 fi

 (echo --------------------------	) >> ./log.txt

 (echo Start Flow			) >> ./log.txt
 (time ./temp/GridFlow/Flow.exe       ) &> ./log/GridFlow/log_step${i}_Flow.txt
 (time ./temp/GridFlow/ViewFlow.exe   ) &> ./log/GridFlow/log_step${i}_ViewFlow.txt
 cp -r ./result/GridFlow/flow ./result/${RESULTDIR}/flow_$i
 cp -r ./data/GridFlow ./result/${RESULTDIR}/flowdata_$i
 (echo Complete Flow			) >> ./log.txt
 (echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt

 (echo --------------------------	) >> ./log.txt

 #MPS_INITIAL
 if [ $i -eq 1 ]; then
  (echo Start MPS initial		) >> ./log.txt
  (time ./temp/MPSini/MPS_ini.exe     ) &> ./log/MPSini/log_IniMPS.txt
  (echo Complete MPS initial		) >> ./log.txt
  (echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt
 fi

 (echo --------------------------	) >> ./log.txt

 #MPS
 (echo Start MPS			) >> ./log.txt
 (time ./temp/MPS/MPS.exe             ) &> ./log/MPS/log_step${i}_MPS.txt
 rsync -r ./result/MPS ./result/${RESULTDIR}
 cp ./result/MPS/restart_file.dat ./result/${RESULTDIR}/MPS/restart_file_${i}.dat
 (echo Complete MPS			) >> ./log.txt
 (echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt

 (echo --------------------------	) >> ./log.txt

 #ICM
 (echo Start ICM			) >> ./log.txt
 (time ./temp/ICM/ICM.exe             ) &> ./log/ICM/log_step${i}_ICM.txt
 cp -r ./result/ICM ./result/${RESULTDIR}/ICM_$i
 (echo Complete MPS			) >> ./log.txt
 (echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt

done

(echo ==========================	) >> ./log.txt
(echo Computation END			) >> ./log.txt
(echo `date '+%Y-%m-%d %H:%M:%S'`	) >> ./log.txt
(echo ==========================	) >> ./log.txt

cp ./log.txt ./result/${RESULTDIR}/log.txt
