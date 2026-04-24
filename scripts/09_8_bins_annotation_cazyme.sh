#!/bin/bash

# set paths
base_dir="path to a working directory"
scriptpath="$base_dir/bins/das_ani98/dbcan-scripts"
dbcanpath="$base_dir/bins/das_ani98/dbcan"

mkdir -p "$scriptpath"
mkdir -p "$dbcanpath"

# Loop over Bakta outputs
for bin in "$base_dir"/bins/das_ani98/bakta/*; do
    b=$(basename "$bin")

    faa="$bin/${b}.faa"
    gff="$bin/${b}.gff3"

    # Skip if files are missing (safety)
    [[ -f "$faa" && -f "$gff" ]] || {
        echo "Skipping $b (missing faa/gff3)"
        continue
    }

    script_name="${scriptpath}/run_${b}.sh"
    log_name="${scriptpath}/${b}.log"
    err_name="${scriptpath}/${b}.err"

    cat <<EOF > "$script_name"
#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N dbcan-${b}
#PBS -e ${err_name}
#PBS -o ${log_name}
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=24:00:00

# Load modules
module load tools
module load signalp/6.0h-slow_sequential
module load dbcan/4.1.4

mkdir -p $dbcanpath/${b}

# Run dbCAN
run_dbcan \
    $faa protein \
    -c $gff \
    --tools all \
    --dia_cpu 40 \
    --hmm_cpu 40 \
    --tf_cpu 40 \
    --stp_cpu 40 \
    --cgc_substrate \
    --db_dir $base_dir/db/dbcan_db_09072024 \
    --out_dir $dbcanpath/${b}
EOF

    chmod +x "$script_name"
    qsub "$script_name"
done
