#!/bin/bash

/lfs/l7/hess/software/python/miniconda/miniconda/condabin/conda activate hap
source /lfs/l7/hess/users/spencers/realdata/root_build/bin/thisroot.sh
source /lfs/l7/hess/users/spencers/realdata/hap/thishess.sh

python /lfs/l7/hess/users/spencers/realdata/process_params.py 127825 /lfs/l7/hess/users/spencers/realdata
