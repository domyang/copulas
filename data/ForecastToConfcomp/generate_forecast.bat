REM !/usr/bin/python
REM  runs ForecastToConfComp

python ForecastToConfComp.py ^
--time-of-forecast=0 ^
--name-of-forecast-file=processed_EMS_forecast_11-15--01-16.csv ^
--positive-boundary=2 ^
--negative-boundary=-2 ^
--critical-number-of-times=4 ^
--intraday-forecast ^
--list-of-input-file-name="list_of_my_input_files.dat" ^
--input-directory="..\..\..\..\..\bpa\from_melanie\input_for_ForecastToConfComp" ^
--output-directory="..\..\..\..\..\bpa\from_melanie\ProcessedForecast"