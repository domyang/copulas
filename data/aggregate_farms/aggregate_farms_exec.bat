REM !/usr/bin/bash

python aggregate_farms.py ^
--data "..\..\..\..\..\bpa\from_melanie\Wind data 1-5-16.xlsx" ^
--vendor EMS ^
--actuals-output-directory "..\..\..\..\..\bpa\from_melanie\actuals_for_prescient" ^
--actuals-writing-file actuals_1-5-16.csv ^
--forecast-output-directory "..\..\..\..\..\bpa\from_melanie\input_for_ForecastToConfComp" ^
--forecast-writing-file EMS_forecast_1-5-16.csv ^
--lower-limit-output-directory "..\..\..\..\..\bpa\from_melanie\input_for_VendorCiToConfComp" ^
--lower-limit-writing-file EMS_lower_limit_1-5-16.csv ^
--upper-limit-output-directory "..\..\..\..\..\bpa\from_melanie\input_for_VendorCiToConfComp" ^
--upper-limit-writing-file EMS_upper_limit_1-5-16.csv


REM --gen-data "..\..\..\..\..\bpa\from_scott\gen_data\2015 09.csv" ^
REM --forecast-mean "..\..\..\..\..\bpa\from_scott\forecast_data\Oct-6-2015\OFFICIAL_gen_mw.txt" ^
REM --lower-limit "..\..\..\..\..\bpa\from_scott\forecast_data\Oct-6-2015\OFFICIAL_gen_min_mw.txt" ^
REM --upper-limit "..\..\..\..\..\bpa\from_scott\forecast_data\Oct-6-2015\OFFICIAL_gen_max_mw.txt" ^
REM --vendor OFFICIAL ^
