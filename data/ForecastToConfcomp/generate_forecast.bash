#!/usr/bin/python
# runs ForecastToConfComp

python ForecastToConfComp.py \
--time-of-forecast=10 \
--name-of-forecast-file=forecast_my.csv \
--positive-boundary=2 \
--negative-boundary=-2 \
--critical-number-of-times=4 \
--intraday-forecast \
--list-of-input-file-name="list_of_my_input_files.dat" \
--input-directory="../../../../../bpa/from_melanie/input_for_ForecastToConfComp" \
--output-directory="../../../../../bpa/from_melanie/ProcessedForecast"