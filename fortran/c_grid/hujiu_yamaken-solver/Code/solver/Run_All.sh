#!/bin/sh
#PBS -l nodes=yato237:ppn=16
#PBS -N ABE_GRIDFLOW_LRE
#PBS -o /dev/null
#PBS -j oe

OMP_NUM_THREADS=16; export OMP_NUM_THREADS
ulimit -s unlimited
cd $PBS_O_WORKDIR

#(echo 6 ; echo 2) | (time ./Grid.exe) &> Log_Grid0.txt
#(time ./IniFlow.exe) &> Log_IniFlow0.txt
#(time ./Flow.exe) &> Log_Flow0.txt
#(time ./ViewFlow.exe) &> Log_ViewFlow0.txt
#(time ./Droplet.exe) &> Log_Droplet0.txt
#(time ./Icing.exe) &> Log_Icing0.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh0.txt

#(time ./Remesh.exe) &> Log_Remesh0.txt
#(time ./Validation.exe) &> Log_Validation0.txt
#(time ./Flow.exe) &> Log_Flow1.txt
#(time ./ViewFlow.exe) &> Log_ViewFlow1.txt
#(time ./Droplet.exe) &> Log_Droplet1.txt
#(time ./Icing.exe) &> Log_Icing1.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh1.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh1.txt
#(time ./Validation.exe) &> Log_Validation1.txt
#(time ./Flow.exe) &> Log_Flow2.txt
#(time ./ViewFlow.exe) &> Log_ViewFlow2.txt
#(time ./Droplet.exe) &> Log_Droplet2.txt
#(time ./Icing.exe) &> Log_Icing2.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh2.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh2.txt
#(time ./Validation.exe) &> Log_Validation2.txt
#(time ./Flow.exe) &> Log_Flow3.txt
#(time ./ViewFlow.exe) &> Log_ViewFlow3.txt
#(time ./Droplet.exe) &> Log_Droplet3.txt
#(time ./Icing.exe) &> Log_Icing3.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh3.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh3.txt
#(time ./Validation.exe) &> Log_Validation3.txt
#(time ./Flow.exe) &> Log_Flow4.txt
#(time ./ViewFlow.exe) &> Log_ViewFlow4.txt
#(time ./Droplet.exe) &> Log_Droplet4.txt
#(time ./Icing.exe) &> Log_Icing4.txt

#(echo 6 ; echo 2) | (time ./Remesh.exe) &> Log_Remesh4.txt
#(time ./Validation.exe) &> Log_Validation4.txt


(echo 6 ; echo 2) | (time ./Grid.exe) &> Log_Grid0.txt
(time ./IniFlow.exe) &> Log_IniFlow0.txt
for i in `seq -w 0 4`
do
(time ./Flow.exe) &> Log_Flow${i}.txt
(time ./ViewFlow.exe) &> Log_ViewFlow${i}.txt
(time ./Droplet.exe) &> Log_Droplet${i}.txt
(time ./Icing.exe) &> Log_Icing${i}.txt
(time ./Remesh.exe) &> Log_Remesh${i}.txt
(time ./Validation.exe) &> Log_Validation${i}.txt
done

