#! /bin/bash
python run.py --seg-N=6 \
                 --seg-s=0.1 \
                 --seg-kappa=50 \
                 --L1Linf-solver=gurobi \
                 --L2Norm-solver=gurobi \
                 --epifit-error-norm=L2 \   

