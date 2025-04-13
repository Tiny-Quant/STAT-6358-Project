# 

## Sampling 

1. Pull the SRA data from the NIH and generate a `.fastq` file for each sample. 

```bash
cd data/
bash pull_SRA_and_fastq.sh 
```

2. Create the necessary alignment indices. 

```bash
cd data/align_indices
bash create align_indices.sh 
```

3. Run the sampling script $n$ times. 

```bash
```

