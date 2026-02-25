# Snakemake Batch Processing Pipeline

A Snakemake-based data pipeline extended to support batch processing of multiple datasets in parallel, designed for execution on a high-performance computing (HPC) cluster.

---

## Overview

This project builds on an existing Snakemake pipeline with the following additions:

- **Batch processing (in testing)** — users can now provide a list of project IDs as input, and the pipeline will process all datasets in a single run
- **HPC parallelization (Coming soon)** — jobs are distributed and executed in parallel across an HPC cluster, significantly reducing processing time at scale
- **Version control** — source code and workflow changes are managed using Git

---

## Requirements

- Python 3.x
- Snakemake
- Access to an HPC cluster
- Git

---

## Usage

### Input

Provide a list of project IDs in a configuration file or as a comma-separated argument:

```yaml
# config.yaml
project_ids:
  - project_001
  - project_002
  - project_003
```

### Running Locally

```bash
snakemake --configfile config.yaml
```

### Running on HPC Cluster

```bash
snakemake --configfile config.yaml --cluster "sbatch --ntasks={threads}" --jobs 10
```

The `--jobs` flag controls how many jobs run in parallel on the cluster.

---

## What Was Added

### Batch Processing (Testing)
Previously, the pipeline accepted a single dataset as input. It has been rewritten to accept a list of project IDs, iterating over each and applying the full pipeline workflow to every dataset automatically.

### HPC Parallel Execution (Coming Soon)
The pipeline is now configured to distribute jobs across an HPC cluster using Snakemake's cluster execution mode. Each project ID is processed as an independent job, allowing multiple datasets to be handled simultaneously.

---

## Version Control

This project is managed with Git. To clone the repository:

```bash
git clone https://github.com/wolfworldrun/<repo-name>
```

---

## Author

Jaliyah Harrison  
[GitHub](https://github.com/wolfworldrun) • jaliyahdharrison@gmail.com
