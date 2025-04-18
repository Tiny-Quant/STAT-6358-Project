# 

## Reproducible Environment

### Single Machine 
1. Ensure docker and docker-compose is installed. 

```bash
docker --version
# Docker version 28.0.4, build b8034c0
docker-compose --version
# Docker Compose version v2.34.0-desktop.1
```

2. Spin up the container. 
```bash 
cd .devcontainer/ 
docker-compose build docker-compose.yml
docker-compose up docker-compose.yml 
docker-compose exec bash 
```
This should drop you into a bash terminal where all the dependencies are installed. 
Ensure that the conda environment `stat6358` is activated. 

### HPC Cluster
Typically HPCs don't allow pure docker containers to run because of security issues.  
Most will allows containers converted to a `singularity` image to run via `apptainer`.

1. Convert the container into a `singularity` image. This has to be run on a machine that has docker installed.

```bash
docker --version
# Docker version 28.0.4, build b8034c0
docker-compose --version
# Docker Compose version v2.34.0-desktop.1
```

```bash
cd .devcontainer/
bash build_image.sh 
```

Now moving to the HPC. 

2. Change the `img_directory` variable in `devcontainer.lua` to the path on the HPC where you have cloned the repo. Make sure it points to the `.devcontainer` folder, not the repo parent directory.  

3. Copy the `.sif` file into `.devcontainer/` on the HPC.  

4. To drop into the container bash shell run: 
```bash
module use /path_on_HPC/STAT-6358-Project/.devcontainer/
module load devcontainer
cd /path_on_HPC/STAT-6358-Project/
singularity shell /path_on_HPC/STAT-6358-Project/.devcontainer/{sif_file_name}.sif
```

5. See `/HPC/template.sbatch` for an example of executing code in the container on an HPC via slurm. 

## Sampling 

1. Pull the SRA data from the NIH and generate a `.fastq` file for each sample. 

```bash
cd .data/
bash pull_SRA_and_fastq.sh 
```

2. Create the necessary alignment indices. 

```bash
cd data/align_indices
bash create align_indices.sh 
```

3. Run the sampling script. 
```bash
cd HPC/
sbatch slurm_array_template.sbatch
```

3a. If you are not using an HPC, then you can run the script directly. 
```bash
cd data/gen_samples
bash generate_sample.sh 
```
Note that each run produces one sample or one line of data, so loop or run the script in parallel $n$ times. 

4. Run the R script to combine the samples into dataframes and write them to `.csv` files.   
```bash
Rscript data/gen_samples/combine_samples.r
```

Now `count_sd_df.csv` and `DE_sd_df.csv` should appear in the `data/gen_samples/` directory. 
