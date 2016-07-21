#! /bin/bash

python -c \
"import code;
import sys;
sys.path.append('/home/ambroiseidoine/UCD/tests');
import error_correlation as ec;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy as sc;
import dateutil.parser as dt;
from dateutil.relativedelta import relativedelta;
code.interact(local=locals())"
