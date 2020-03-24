#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
MOLOC=$HOME/bin/eczema_gwas_fu/colocalization/moloc
MOLOC_ANALYSIS=$HOME/analysis/colocalization/moloc
#Submit test Moloc run
cd $MOLOC_ANALYSIS
qsub $MOLOC/sub_moloc_test.sh
