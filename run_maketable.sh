#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4g
#SBATCH --output=genomestats.out
perl scripts/make_summary_table.pl > fungi_genome_stats.tab 2> err &
