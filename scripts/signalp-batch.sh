#!/bin/bash

cd /home/projects/dtu_00032/analysis/milk-cohort
scriptpath="bins/new-ani98/signalp-scripts"
SPcanpath="bins/new-ani98/bakta-signalp"
mkdir -p $scriptpath
mkdir -p $SPcanpath


for bin in bins/new-ani98/bakta/*; do
    b=$(basename $bin)
    script_name="${scriptpath}/run_${b}.sh"
    log_name="${scriptpath}/${b}.log"
    err_name="${scriptpath}/${b}.err"

    echo "#!/bin/bash" > $script_name
    echo "#PBS -W group_list=dtu_00032 -A dtu_00032" >> $script_name
    echo "#PBS -N signalp-${b}" >> $script_name
    echo "#PBS -e signalp-${b}.err" >> $script_name
    echo "#PBS -o signalp-${b}.log" >> $script_name
    echo "#PBS -l nodes=1:thinnode:ppn=40" >> $script_name
    echo "#PBS -l mem=188gb" >> $script_name
    echo "#PBS -l walltime=4:00:00" >> $script_name
    echo "" >> $script_name
    echo "module load tools" >> $script_name
    echo "module load signalp/6.0h-slow_sequential  " >> $script_name
    echo "cd /home/projects/dtu_00032/analysis/milk-cohort" >> $script_name
    echo "" >> $script_name
    echo "mkdir -p $SPcanpath/$b" >> $script_name
    echo "signalp6 \\" >> $script_name
    echo "    --fastafile $bin/$b.faa \\" >> $script_name
    echo "    --output_dir $SPcanpath/$b \\" >> $script_name
    echo "    --organism other \\" >> $script_name
    echo "    --mode slow-sequential \\" >> $script_name
    echo "    -wp 40 ;" >> $script_name

    chmod +x $script_name
    qsub $script_name
done
