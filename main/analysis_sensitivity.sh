#!/bin/bash
#$ -N network_sosd_sensitivity
#$ -j y
#$ -o ../../../out/networkSoSD/qsub/
#$ -l m_mem_free=50G

Rscript --no-save analysis_sensitivity.R
