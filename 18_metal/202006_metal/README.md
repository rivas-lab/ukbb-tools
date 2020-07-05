# UKB meta-analysis

https://github.com/rivas-lab/ukbb-tools/issues/22


```
sbatch -p mrivas --qos=high_p --nodes=1 --mem=4000 --cores=1 --time=3:00:00 --job-name=metal --output=logs/metal.%A_%a.out --error=logs/metal.%A_%a.err --array=1-402 1_metal.sbatch.sh
Submitted batch job 3584555
```