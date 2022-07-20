CONFIG=$1

snakemake -s /project/holstegelab/Software/snakemake_pipeline/bin/Snakefile_merge_and_pileup.py --latency-wait 60 --configfile ${CONFIG} --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config /project/holstegelab/Software/snakemake_pipeline/config/config_cluster_merge_and_pileup.yml --jobs 1
