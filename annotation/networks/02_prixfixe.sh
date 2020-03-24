#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/networks/prixfixe
scripts=$HOME/bin/eczema_gwas_fu/annotation/networks

cd $analysis

Rscript $scripts/prixfixe.R