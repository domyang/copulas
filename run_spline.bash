#! /bin/bash
python spline.py --seg-N=5 \
                 --seg-s=0.1 \
                 --seg-kappa=50 \
                 --L1Linf-solver=gurobi \
                 --L2Norm-solver=gurobi \
                 --epifit-error-norm=L2 \   

