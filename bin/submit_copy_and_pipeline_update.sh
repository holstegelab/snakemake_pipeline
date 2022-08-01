CONFIG=$1

snakemake -s /project/holstegelab/Software/snakemake_pipeline/bin/Snakefile_with_copy_update.py --latency-wait 60 --configfile ${CONFIG} --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config /project/holstegelab/Software/snakemake_pipeline/config/config_cluster.yml --jobs 1