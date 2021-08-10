#!/bin/bash
#$ -N simul_switching
#$ -j y
#$ -o ../results/qsub/
#$ -l m_mem_free=12G
#$ -pe openmpi 4

Rscript --no-save simulation_rho_switching.R
