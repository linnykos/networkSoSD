#!/bin/bash
#$ -N network_sosd_sensitivity
#$ -j y
#$ -o .../../../out/networkSoSD/qsub/
#$ -l m_mem_free=15G

Rscript --no-save analysis_sensitivity.R
