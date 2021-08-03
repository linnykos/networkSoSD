#!/bin/bash
#$ -N simulation_rho
#$ -j y
#$ -o ../results/qsub/
#$ -l m_mem_free=12G
#$ -pe openmpi 4

Rscript --no-save simulation_rho.R
