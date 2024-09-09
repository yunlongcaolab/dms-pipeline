#!/bin/bash

set -e

if [ -z $1 ]; then
    echo "Please provide the project name as input!"
    exit 1
fi

workdir=/lustre/grp/cyllab/share/DMS/$1

python ext_utils/merge_escape_scores.py -i $workdir/processed -o $workdir/merge