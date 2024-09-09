#!/bin/bash

workdir=/lustre/grp/cyllab/share/DMS/$1

python ext_utils/merge_escape_scores.py -i $workdir/processed -o $workdir/merge