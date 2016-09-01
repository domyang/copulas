#!/usr/bin/python

python VendorCiToConfComp.py \
--time-of-forecast=0 \
--name-of-lower-limit-file=processed_EMS_lower_limit_11-15--01-16.csv \
--name-of-upper-limit-file=processed_EMS_upper_limit_11-15--01-16.csv \
--intraday-forecast \
--list-of-input-file-name="list_of_my_input_files.dat" \
--input-directory="../../../../../bpa/from_melanie/input_for_VendorCiToConfComp" \
--lower-limit-output-directory="../../../../../bpa/from_melanie/ProcessedLowerLimit" \
--upper-limit-output-directory="../../../../../bpa/from_melanie/ProcessedUpperLimit"
