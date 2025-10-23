# this instructs snakemake to submit each rule as a separate job with specified default and job specific resources (which are defined in the cluster.yaml file) 

snakemake   \
        -s ../Snakefile --jobs 100  --use-conda --conda-frontend conda  --cluster \
        "sbatch -p {cluster.partition} --job-name={rule} --mem={cluster.mem} --time={cluster.time} --ntasks={cluster.ntasks} --nodes={cluster.nodes}  --output=slurm-logs/%x.%j.log --verbose" \
        --cluster-config ../cluster.yaml --latency-wait 120
