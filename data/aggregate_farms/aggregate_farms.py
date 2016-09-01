from optparse import OptionParser
import openpyxl
import os
import math
import csv
from datetime import datetime, timedelta
try:
    from collections import OrderedDict
except:
    from ordereddict import OrderedDict
from ucd_to_pst import get_daylight_saving_time, get_daylight_saving_end


parser = OptionParser()
parser.add_option('--list-of-data-files', help='dat-file with a list of the xlsx-files with the forecasts and actuals', action='store', dest='list_of_data_files', type='string', default=None)
parser.add_option('--input-directory', dest='input_directory', action='store', type='string', default='input_files')
parser.add_option('--gen-data', help='csv-file with the actuals', action='store', dest='gen_data', type='string', default=None)
parser.add_option('--forecast-mean', help='txt-file with the forecasted means ', action='store', dest='forecast_mean', type='string', default=None)
parser.add_option('--lower-limit', help='txt-file with the lower limits', action='store', dest='lower_limit', type='string', default=None)
parser.add_option('--upper-limit', help='txt-file with the upper limits', action='store', dest='upper_limit', type='string', default=None)
parser.add_option('--forecast-writing-file', help='csv-file where the forecast has to be stored', action='store', dest='forecast_writing_file', type='string', default=None)
parser.add_option("--forecast-output-directory", help="The directory where the forecast output data is stored.The default is: 'output_files'", action="store", dest="forecast_output_directory", type="string", default="output_files")
parser.add_option('--lower-limit-writing-file', help='csv-file where the lower limits has to be stored', action='store', dest='lower_limit_writing_file', type='string', default=None)
parser.add_option("--lower-limit-output-directory", help="The directory where the lower limit output data is stored.The default is: 'output_files'", action="store", dest="lower_limit_output_directory", type="string", default="output_files")
parser.add_option('--upper-limit-writing-file', help='csv-file where the upper limits has to be stored', action='store', dest='upper_limit_writing_file', type='string', default=None)
parser.add_option("--upper-limit-output-directory", help="The directory where the upper limit output data is stored.The default is: 'output_files'", action="store", dest="upper_limit_output_directory", type="string", default="output_files")
parser.add_option('--actuals-writing-file', help='csv-file where the actuals has to be stored', action='store', dest='actuals_writing_file', type='string', default=None)
parser.add_option("--actuals-output-directory", help="The directory where the actual output data is stored.The default is: 'output_files'", action="store", dest="actuals_output_directory", type="string", default="output_files")
parser.add_option('--vendor', help='vendor', action='store', dest='vendor', type='string', default='OFFICIAL')
parser.add_option('--max-nb-interpolated-hours', dest='max_nb_interpolated_hours', type='float', default=2, help='maximal numbers of hours in a gap that should be interpolated')
(options, args) = parser.parse_args()

if options.max_nb_interpolated_hours < 0 or options.max_nb_interpolated_hours > 23:
    raise RuntimeError('max-nb-interpolated-hours has to be an interger greater than 0 and smaller than 24!')

days = []
times = []
days_utc = []
times_utc = []

if options.list_of_data_files != None:
    #read file with files
    input_file_list = []
    list_file = open(options.list_of_data_files, 'r')
    for line in list_file:
        if line[0] == '#':
            continue
        name_of_file = str(line.strip('\n'))
        if name_of_file in input_file_list:
            continue
        input_file_list.append(name_of_file)
    list_file.close()

    actuals = OrderedDict()     # {'start time of observation (1h)': {'farm': 'observation'}}
    farms_actual = []
    forecast = OrderedDict()    # {datetime of forecast: {Hr: {farm: forecast}}}
    farms_forecast = []
    lower_limits = OrderedDict()    # {datetime of forecast: {Hr: {farm: lower_limit}}}
    farms_lower_limit = []
    upper_limits = OrderedDict()    # {datetime of forecast: {Hr: {farm: upper_limit}}}         
    farms_upper_limit = []    

    for file_name in input_file_list:
        name_of_file = options.input_directory + os.sep + file_name
        data = openpyxl.load_workbook(name_of_file, 'r')
        ########################
        ##### read actuals #####
        ########################
        print('Reading actuals..')
        if 'Actual' not in data.get_sheet_names():
            raise RuntimeError('ERROR:', file_name, 'has no sheet named Actuals!')
        actuals_sheet = data['Actual']
        IndexFarmMap = OrderedDict()
        for line in actuals_sheet.rows:
            line = line[0].value
            if line == None:
                continue
            actual_line = line.strip('\n').split(',')
            if actual_line[0] == 'Date/Time (UTC)':
                for i in range(2, len(actual_line)):
                    IndexFarmMap[i] = actual_line[i].strip()
                    if actual_line[i].strip() not in farms_actual:
                        farms_actual.append(actual_line[i].strip())
            else:
                try:
                    start_time_utc = datetime.strptime(actual_line[0], '%Y-%m-%d %H:%M:%S')
                    if start_time_utc.date() not in days_utc:
                        days_utc.append(start_time_utc.date())
                    if start_time_utc.time() not in times_utc:
                        times_utc.append(start_time_utc.time())
                    daylight_saving_start = get_daylight_saving_time(str(start_time_utc.year))
                    daylight_saving_end = get_daylight_saving_end(str(start_time_utc.year))
                    difference_pdt_tu_utc = 'will be set now'
                    if start_time_utc >= daylight_saving_start and start_time_utc <= daylight_saving_end:
                        difference_pdt_to_utc = 7
                    elif start_time_utc < daylight_saving_start or start_time_utc > daylight_saving_end:
                        difference_pdt_to_utc = 8
                    start_time = start_time_utc - timedelta(hours=difference_pdt_to_utc)
                    if start_time.date() not in days:
                        days.append(start_time.date())
                    if start_time.time() not in times:
                        times.append(start_time.time())
                    if start_time not in actuals.keys():
                        actuals[start_time] = OrderedDict()
                    for i in range(2, len(actual_line)):
                        farm = IndexFarmMap[i]
                        actuals[start_time][farm] = actual_line[i].strip()
                except ValueError:
                    continue
    
        Vendor = options.vendor
        if Vendor == 'OFFICIAL':
            vendor = 'Official'
        else:
            vendor = Vendor
        
        #########################
        ##### read forecast #####
        #########################
        print('Reading forecasts..')
        if vendor + ' Mean' not in data.get_sheet_names():
            raise RuntimeError('ERROR:', options.data_file, 'has no sheet named', vendor, 'Mean', '!')
        mean_sheet = data[vendor + ' Mean']
        for line in mean_sheet.rows:
            line = line[0].value
            if line == None:
                continue
            actual_line = line.strip('\n').split(',')
            if actual_line[0] == Vendor:
                dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
                farm = actual_line[2]
                if farm not in farms_forecast:
                    farms_forecast.append(farm)
                if dt_utc not in forecast.keys():
                    forecast[dt_utc] = OrderedDict()
                for i in range(168):
                    if i+1 not in forecast[dt_utc].keys():
                        forecast[dt_utc][i+1] = OrderedDict()
                    mean = actual_line[i+3]
                    forecast[dt_utc][i+1][farm] = mean
                
        #############################
        ##### read lower limits #####
        #############################
        print('Reading lower limits..')
        if vendor + ' Min' not in data.get_sheet_names():
            raise RuntimeError('ERROR:', options.data_file,' has no sheet named', vendor, 'Min', '!')
        min_sheet = data[vendor + ' Min']
        for line in min_sheet.rows:
            line = line[0].value
            if line == None:
                continue
            actual_line = line.strip('\n').split(',')
            if actual_line[0] == Vendor:
                dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
                farm = actual_line[2]
                if farm not in farms_lower_limit:
                    farms_lower_limit.append(farm)
                if dt_utc not in lower_limits.keys():
                    lower_limits[dt_utc] = OrderedDict()
                for i in range(168):
                    if i+1 not in lower_limits[dt_utc].keys():
                        lower_limits[dt_utc][i+1] = OrderedDict()
                    lower_limit = actual_line[i+3]
                    lower_limits[dt_utc][i+1][farm] = lower_limit
                    
        #############################
        ##### read upper limits #####
        #############################
        print('Reading upper limits..')
        if vendor + ' Max' not in data.get_sheet_names():
            raise RuntimeError('ERROR: ', options.data_file, ' has no sheet named', vendor, 'Max', '!')
        max_sheet = data[vendor + ' Max']
        for line in max_sheet.rows:
            line = line[0].value
            if line == None:
                continue
            actual_line = line.strip('\n').split(',')
            if actual_line[0] == Vendor:
                dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
                farm = actual_line[2]
                if farm not in farms_upper_limit:
                    farms_upper_limit.append(farm)
                if dt_utc not in upper_limits.keys():
                    upper_limits[dt_utc] = OrderedDict()
                for i in range(168):
                    if i+1 not in upper_limits[dt_utc].keys():
                        upper_limits[dt_utc][i+1] = OrderedDict()
                    upper_limit = actual_line[i+3]
                    upper_limits[dt_utc][i+1][farm] = upper_limit

elif options.gen_data!=None and options.forecast_mean!=None and options.lower_limit!=None and options.upper_limit!=None:
    ########################
    ##### read actuals #####
    ########################
    print('Reading actuals..')
    with open(options.gen_data) as gen_data:
        actuals_reader = csv.reader(gen_data, delimiter=',', quotechar = '|')
        actuals = OrderedDict()     # {'start time of observation (1h)': {'farm': 'observation'}}
        farms_actual=[]
        IndexFarmMap = OrderedDict()
        for line in actuals_reader:
            if line[0] == '' or line[0] == ' ':
                continue
            elif line[0] == 'start time':
                for i in range(1, len(line)):
                    IndexFarmMap[i] = line[i]
                    if line[i].strip() not in farms_actual:
                        farms_actual.append(line[i].strip())
            else:
                start_time = datetime.strptime(line[0], '%d-%b-%y %H:%M:%S')
                if start_time.date() not in days:
                    days.append(start_time.date())
                if start_time.time() not in times:
                    times.append(start_time.time())
                daylight_saving_start = get_daylight_saving_time(str(start_time.year))
                daylight_saving_end = get_daylight_saving_end(str(start_time.year))
                difference_pdt_tu_utc = 'will be set now'
                if start_time >= daylight_saving_start and start_time <= daylight_saving_end:
                    difference_pdt_to_utc = 7
                elif start_time < daylight_saving_start or start_time > daylight_saving_end:
                    difference_pdt_to = 8
                start_time_utc = start_time + timedelta(hours=difference_pdt_to_utc)
                if start_time_utc.date() not in days_utc:
                    days_utc.append(start_time_utc.date())
                if start_time_utc.time() not in times_utc:
                    times_utc.append(start_time_utc.time())
                if start_time not in actuals.keys():
                    actuals[start_time] = OrderedDict()
                for i in range(1, len(line)):
                    farm = IndexFarmMap[i]
                    actuals[start_time][farm] = line[i]
    
    #########################
    ##### read forecast #####
    #########################
    print('Reading forecasts..')
    mean_file = open(options.forecast_mean, 'r')
    forecast = OrderedDict()    # {datetime of forecast: {Hr: {farm: forecast}}}
    farms_forecast=[]
    for line in mean_file:
        actual_line = line.strip('\n').split(',')
        if actual_line[0] == options.vendor:
            dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
            farm = actual_line[2]
            if farm not in farms_forecast:
                farms_forecast.append(farm)
            if dt_utc not in forecast.keys():
                forecast[dt_utc] = OrderedDict()
            for i in range(168):
                if i+1 not in forecast[dt_utc].keys():
                    forecast[dt_utc][i+1] = OrderedDict()
                mean = actual_line[i+3]
                forecast[dt_utc][i+1][farm] = mean
                
    #############################
    ##### read lower_limits #####
    #############################
    print('Reading lower limits..')
    min_file = open(options.lower_limit, 'r')
    lower_limits = OrderedDict()    # {datetime of forecast: {Hr: {farm: lower_limit}}}
    farms_lower_limit=[]
    for line in min_file:
        actual_line = line.strip('\n').split(',')
        if actual_line[0] == options.vendor:
            dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
            farm = actual_line[2]
            if farm not in farms_lower_limit:
                farms_lower_limit.append(farm)
            if dt_utc not in lower_limits.keys():
                lower_limits[dt_utc] = OrderedDict()
            for i in range(168):
                if i+1 not in lower_limits[dt_utc].keys():
                    lower_limits[dt_utc][i+1] = OrderedDict()
                min = actual_line[i+3]
                lower_limits[dt_utc][i+1][farm] = min
                
    #############################
    ##### read upper_limits #####
    #############################
    print('Reading upper limits..')
    max_file = open(options.upper_limit, 'r')
    upper_limits = OrderedDict()    # {datetime of forecast: {Hr: {farm: upper_limit}}}
    farms_upper_limit=[]
    for line in max_file:
        actual_line = line.strip('\n').split(',')
        if actual_line[0] == options.vendor:
            dt_utc = datetime.strptime(actual_line[1], '%Y-%m-%d %H:%M:%S')
            farm = actual_line[2]
            if farm not in farms_upper_limit:
                farms_upper_limit.append(farm)
            if dt_utc not in upper_limits.keys():
                upper_limits[dt_utc] = OrderedDict()
            for i in range(168):
                if i+1 not in upper_limits[dt_utc].keys():
                    upper_limits[dt_utc][i+1] = OrderedDict()
                max = actual_line[i+3]
                upper_limits[dt_utc][i+1][farm] = max

else:
    raise RuntimeError('Not enough input data!')


#############################
##### find common farms #####
#############################
print('Finding common farms..')
list_of_farms = []
for farm in farms_actual:
    flag = 1
    if farm not in farms_forecast:
        print('There are no forecasts for farm ' + farm + '!')
        flag = 0
    if farm not in farms_lower_limit:
        print('There are no lower limits for farm ' + farm + '!')
        flag = 0
    if farm not in farms_upper_limit:
        print('There are no upper limits for farm ' + farm + '!')
        flag = 0
    if flag == 1:
        list_of_farms.append(farm)
string_of_farms = ''
for farm in list_of_farms:
    string_of_farms = string_of_farms + farm + ', '
print('Considered farms:', string_of_farms.strip(', '))


#####################################
##### interpolate missing hours #####
#####################################
#actuals
print('Interpolating actuals..')
for d in days:
    for t in times:
        dt = datetime.combine(d, t)
        if dt not in actuals.keys():
            print('***WARNING: Actual is missing for date', dt.date(), 'and time', dt.time(), 'for all farms. Interpolating using one hour increments...')
            actuals[dt] = OrderedDict()
            flag = 1
            failed_interpolation = ''
            for farm in list_of_farms:
                dti, i = dt, 1
                while True: # We find the closest lower neighbor
                    dti = dt - timedelta(minutes=i*60)
                    if dti in actuals.keys() and farm in actuals[dti].keys() and actuals[dti][farm]!='':
                        break
                    else:
                        i = i+1
                    if i > options.max_nb_interpolated_hours : # if we don't have data for 24 hours there is a bigger issue
                        dti = -1
                        break
                dtf, f = dt, 1
                while True: # We find the closest upper neighbor
                    dtf = dt + timedelta(minutes=f*60)
                    if dtf in actuals.keys() and farm in actuals[dtf].keys() and actuals[dtf][farm]!='':
                        break
                    else:
                        f = f+1
                    if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                        dtf = -1
                        break
                if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                    flag = 0
                    failed_interpolation = failed_interpolation + farm + ', '
                    continue
                else:
                    actuals[dt][farm] = 0.5*(float(actuals[dti][farm]) + float(actuals[dtf][farm]))
            if flag == 0:
                print("***WARNING: Interpolation failed - we can't find either lower or upper neighbors for date time " + str(dt) + ', within ' + str(options.max_nb_interpolated_hours) + ' hours, for farms: ' + failed_interpolation.strip(', '))
        elif dt in actuals.keys():
            flag1 = 1
            flag2 = 1
            missing_farms = ''
            failed_interpolation = ''
            for farm in list_of_farms:
                if farm not in actuals[dt].keys() or actuals[dt][farm]=='':
                    flag1 = 0
                    missing_farms = missing_farms + farm + ', '
                    dti, i = dt, 1
                    while True: # We find the closest lower neighbor
                        dti = dt - timedelta(minutes=i*60)
                        if dti in actuals.keys() and farm in actuals[dti].keys() and actuals[dti][farm]!='':
                            break
                        else:
                            i = i+1
                        if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dti = -1
                            break
                    dtf, f = dt, 1
                    while True: # We find the closest upper neighbor
                        dtf = dt + timedelta(minutes=f*60)
                        if dtf in actuals.keys() and farm in actuals[dtf].keys() and actuals[dtf][farm]!='':
                            break
                        else:
                            f = f+1
                        if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dtf = -1
                            break
                    if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                        flag2 = 0
                        failed_interpoaltion = failed_interpolation + farm + ', '
                        continue
                    else:
                        actuals[dt][farm] = 0.5*(float(actuals[dti][farm]) + float(actuals[dtf][farm]))
            if flag1 == 0:
                print('***WARNING: Actual is missing fore date ' + str(dt.date()) + ' and time ' + str(dt.time()) + ' for farms: ', missing_farms.strip(', ') + '. Interpolating using one hour increments...')
            if flag2 == 0:
                print("***WARNING: Interpolation failed - we can't find either lower or upper neighbors for date time " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for farms: " + failed_interpolation.strip(', '))

for dt in actuals.keys():
    print(dt)
    for key in actuals[dt].keys():
        print(key)
                
print('')
# forecast
print('Interpolating forecasts..')
for d in days_utc:
    for t in times_utc:
        dt = datetime.combine(d, t)
        if dt not in forecast.keys():            
            print('***WARNING: Forecast is missing for date', dt.date(), 'and time (UTC)', dt.time(), 'for all forecast horizons. Interpolating using one hour increments...')
            forecast[dt] = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                forecast[dt][Hr] = OrderedDict()
                for farm in list_of_farms:
                    dti, i = dt, 1
                    while True: # We find the closest lower neighbor
                        dti = dt - timedelta(minutes=i*60)
                        if dti in forecast.keys() and Hr in forecast[dti].keys() and farm in forecast[dti][Hr].keys() and forecast[dti][Hr][farm]!='':
                            break
                        else:
                            i = i+1
                        if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dti = -1 
                            break
                    dtf, f = dt, 1
                    while True: # We find the closest upper neighbor
                        dtf = dt + timedelta(minutes=f*60)
                        if dtf in forecast.keys() and Hr in forecast[dtf].keys() and farm in forecast[dtf][Hr].keys() and forecast[dtf][Hr][farm]!='':
                            break
                        else:
                            f = f+1
                        if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dtf = -1
                            break
                    if (dti == -1 or dtf == -1) or (i+f-1 > options.max_nb_interpolated_hours):
                        if Hr not in failed_interpolation.keys():
                            failed_interpolation[Hr] = ''
                        failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                        continue
                    else:
                        forecast[dt][Hr][farm] = 0.5*(float(forecast[dti][Hr][farm]) + float(forecast[dtf][Hr][farm]))
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find either lower or upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons
        elif dt in forecast.keys():
            missing_forecast = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                if Hr not in forecast[dt].keys():
                    forecast[dt][Hr] = OrderedDict()
                    missing_forecast[Hr] = ''
                    for farm in list_of_farms:
                        missing_forecast[Hr] = missing_forecast[Hr] + farm + ', '
                        dti, i = dt, 1
                        while True: # We find the closest lower neighbor
                            dti = dt - timedelta(minutes=i*60)
                            if dti in forecast.keys() and Hr in forecast[dti].keys() and farm in forecast[dti][Hr].keys() and forecast[dti][Hr][farm]!='':
                                break
                            else:
                                i = i+1
                            if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dti = -1
                                break
                        dtf, f = dt, 1
                        while True: # We find the closest upper neighbor
                            dtf = dt + timedelta(minutes=f*60)
                            if dtf in forecast.keys() and Hr in forecast[dtf].keys() and farm in forecast[dtf][Hr].keys() and forecast[dtf][Hr][farm]!='':
                                break
                            else:
                                f = f+1
                            if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dtf = -1
                                break
                        if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                            if Hr not in failed_interpolation.keys():
                                failed_interpolation[Hr] = ''
                            failed_interpolation[Hr] = failed_interpoaltion[Hr] + farm + ', '
                            continue
                        else:
                            forecast[dt][Hr][farm] = 0.5*(float(forecast[dti][Hr][farm]) + float(forecast[dtf][Hr][farm]))
                elif Hr in forecast[dt].keys():
                    for farm in list_of_farms:
                        if farm not in forecast[dt][Hr].keys() or forecast[dt][Hr][farm] == '':
                            if Hr not in missing_forecast.keys():
                                missing_forecast[Hr]= ''
                            missing_forecast[Hr] = missing_forecast[Hr] + farm + ', '
                            dti, i = dt, 1
                            while True: # We find the closest lower neighbor
                                dti = dt - timedelta(minutes=i*60)
                                if dti in forecast.keys() and Hr in forecast[dti].keys() and farm in forecast[dti][Hr].keys() and forecast[dti][Hr][farm]!='':
                                    break
                                else:
                                    i = i+1
                                if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dti = -1
                                    break
                            dtf, f = dt, 1
                            while True: # We find the closest upper neighbor
                                dtf = dt + timedelta(minutes=f*60)
                                if dtf in forecast.keys() and Hr in forecast[dtf].keys() and farm in forecast[dtf][Hr].keys() and forecast[dtf][Hr][farm]!='':
                                    break
                                else:
                                    f = f+1
                                if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dtf = -1
                                    break
                            if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                                if Hr not in failed_interpolation.keys():
                                    failed_interpolation[Hr] = ''
                                failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                                continue
                            else:
                                forecast[dt][Hr][farm] = 0.5*(float(forecast[dti][Hr][farm]) + float(forecast[dtf][Hr][farm]))
            visited_hr = []
            for Hr1 in missing_forecast.keys():
                if Hr1 in visited_hr:
                    continue
                horizons = [Hr1]
                for Hr2 in missing_forecast.keys():
                    if Hr2 != Hr1 and missing_forecast[Hr2] == missing_forecast[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print('***WARNING: Forecast is missing for date ' + str(dt.date()) +  ' and time (UTC) ' + str(dt.time()) + ' for forecast horizon(s) ' +  hor.strip(', ') +  ' for farm(s) ' + missing_forecast[Hr1].strip(', ') + '. Interpolating using one hour increments...')
                visited_hr = visited_hr + horizons
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find either lower or upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons

print('')
# lower limits
print('Interpolating the lower limits..')
for d in days_utc:
    for t in times_utc:
        dt = datetime.combine(d, t)
        if dt not in lower_limits.keys():
            print('***WARNING: Lower limit is missing for date', dt.date(), 'and time (UTC)', dt.time(), 'for all forecast horizons. Interpolating using one hour increments...')
            lower_limits[dt] = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                lower_limits[dt][Hr] = OrderedDict()
                for farm in list_of_farms:
                    dti, i = dt, 1
                    while True: # We find the closest lower neighbor
                        dti = dt - timedelta(minutes=i*60)
                        if dti in lower_limits.keys() and Hr in lower_limits[dti].keys() and farm in lower_limits[dti][Hr].keys() and lower_limits[dti][Hr][farm]!='':
                            break
                        else:
                            i = i+1
                        if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dti = -1
                            break
                    dtf, f = dt, 1
                    while True: # We find the closest upper neighbor
                        dtf = dt + timedelta(minutes=f*60)
                        if dtf in lower_limits.keys() and Hr in lower_limits[dtf].keys() and farm in lower_limits[dtf][Hr].keys() and lower_limits[dtf][Hr][farm]!='':
                            break
                        else:
                            f = f+1
                        if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dtf = -1
                            break
                    if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                        if Hr not in failed_interpolation.keys():
                            failed_interpolation[Hr] = ''
                        failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                        continue
                    else:
                        lower_limits[dt][Hr][farm] = 0.5*(float(lower_limits[dti][Hr][farm]) + float(lower_limits[dtf][Hr][farm]))
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find neither lower nor upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons
        elif dt in lower_limits.keys():
            missing_lower_limits = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                if Hr not in lower_limits[dt].keys():
                    lower_limits[dt][Hr] = OrderedDict()
                    missing_lower_limits[Hr] = ''
                    for farm in list_of_farms:
                        missing_lower_limits[Hr] = missing_lower_limits[Hr] + farm + ', '
                        dti, i = dt, 1
                        while True: # We find the closest lower neighbor
                            dti = dt - timedelta(minutes=i*60)
                            if dti in lower_limits.keys() and Hr in lower_limits[dti].keys() and farm in lower_limits[dti][Hr].keys() and lower_limits[dti][Hr][farm]!='':
                                break
                            else:
                                i = i+1
                            if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dti = -1
                                break
                        dtf, f = dt, 1
                        while True: # We find the closest upper neighbor
                            dtf = dt + timedelta(minutes=f*60)
                            if dtf in lower_limits.keys() and Hr in lower_limits[dtf].keys() and farm in lower_limits[dtf][Hr].keys() and lower_limits[dtf][Hr][farm]!='':
                                break
                            else:
                                f = f+1
                            if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dtf = -1
                                break
                        if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                            if Hr not in failed_interpolation.keys():
                                failed_interpolation[Hr] = ''
                            failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                            continue
                        else:
                            lower_limits[dt][Hr][farm] = 0.5*(float(lower_limits[dti][Hr][farm]) + float(lower_limits[dtf][Hr][farm]))
                elif Hr in lower_limits[dt].keys():
                    for farm in list_of_farms:
                        if farm not in lower_limits[dt][Hr].keys() or lower_limits[dt][Hr][farm] == '':
                            if Hr not in missing_lower_limits.keys():
                                missing_lower_limits[Hr]= ''
                            missing_lower_limits[Hr] = missing_lower_limits[Hr] + farm + ', '
                            dti, i = dt, 1
                            while True: # We find the closest lower neighbor
                                dti = dt - timedelta(minutes=i*60)
                                if dti in lower_limits.keys() and Hr in lower_limits[dti].keys() and farm in lower_limits[dti][Hr].keys() and lower_limits[dti][Hr][farm]!='':
                                    break
                                else:
                                    i = i+1
                                if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dti = -1
                                    break
                            dtf, f = dt, 1
                            while True: # We find the closest upper neighbor
                                dtf = dt + timedelta(minutes=f*60)
                                if dtf in lower_limits.keys() and Hr in lower_limits[dtf].keys() and farm in lower_limits[dtf][Hr].keys() and lower_limits[dtf][Hr][farm]!='':
                                    break
                                else:
                                    f = f+1
                                if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dtf = -1
                                    break
                            if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                                if Hr not in failed_interpolation.keys():
                                    failed_interpolation[Hr] = ''
                                failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                                continue
                            else:
                                lower_limits[dt][Hr][farm] = 0.5*(float(lower_limits[dti][Hr][farm]) + float(lower_limits[dtf][Hr][farm]))
            visited_hr = []
            for Hr1 in missing_lower_limits.keys():
                if Hr1 in visited_hr:
                    continue
                horizons = [Hr1]
                for Hr2 in missing_lower_limits.keys():
                    if Hr2 != Hr1 and missing_lower_limits[Hr2] == missing_lower_limits[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print('***WARNING: Lower limit is missing for date ' + str(dt.date()) +  ' and time (UTC) ' + str(dt.time()) + ' for forecast horizon(s) ' +  hor.strip(', ') +  ' for farm(s) ' + missing_lower_limits[Hr1].strip(', ') + '. Interpolating using one hour increments...')
                visited_hr = visited_hr + horizons
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find either lower or upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons

print('')
# upper limits
print('Interpolating the upper limits..')
for d in days_utc:
    for t in times_utc:
        dt = datetime.combine(d, t)
        if dt not in upper_limits.keys():
            print('***WARNING: Upper limit is missing for date', dt.date(), 'and time (UTC)', dt.time(), 'for all forecast horizons. Interpolating using one hour increments...')
            upper_limits[dt] = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                upper_limits[dt][Hr] = OrderedDict()
                for farm in list_of_farms:
                    dti, i = dt, 1
                    while True: # We find the closest lower neighbor
                        dti = dt - timedelta(minutes=i*60)
                        if dti in upper_limits.keys() and Hr in upper_limits[dti].keys() and farm in upper_limits[dti][Hr].keys() and upper_limits[dti][Hr][farm]!='':
                            break
                        else:
                            i = i+1
                        if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dti = -1
                            break
                    dtf, f = dt, 1
                    while True: # We find the closest upper neighbor
                        dtf = dt + timedelta(minutes=f*60)
                        if dtf in upper_limits.keys() and Hr in upper_limits[dtf].keys() and farm in upper_limits[dtf][Hr].keys() and upper_limits[dtf][Hr][farm]!='':
                            break
                        else:
                            f = f+1
                        if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                            dtf = -1
                            break
                    if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                        if Hr not in failed_interpolation.keys():
                            failed_interpolation[Hr] = ''
                        failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                        continue
                    else:
                        upper_limits[dt][Hr][farm] = 0.5*(float(upper_limits[dti][Hr][farm]) + float(upper_limits[dtf][Hr][farm]))
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find neither lower nor upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) +  "hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons
        elif dt in upper_limits.keys():
            missing_upper_limits = OrderedDict()
            failed_interpolation = OrderedDict()
            for Hr in range(1, 169):
                if Hr not in upper_limits[dt].keys():
                    upper_limits[dt][Hr] = OrderedDict()
                    missing_lower_limits[Hr] = ''
                    for farm in list_of_farms:
                        missing_lower_limits[Hr] = missing_upper_limits[Hr] + farm + ', '
                        dti, i = dt, 1
                        while True: # We find the closest lower neighbor
                            dti = dt - timedelta(minutes=i*60)
                            if dti in upper_limits.keys() and Hr in upper_limits[dti].keys() and farm in upper_limits[dti][Hr].keys() and upper_limits[dti][Hr][farm]!='':
                                break
                            else:
                                i = i+1
                            if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dti = -1
                                break
                        dtf, f = dt, 1
                        while True: # We find the closest upper neighbor
                            dtf = dt + timedelta(minutes=f*60)
                            if dtf in upper_limits.keys() and Hr in upper_limits[dtf].keys() and farm in upper_limits[dtf][Hr].keys() and upper_limits[dtf][Hr][farm]!='':
                                break
                            else:
                                f = f+1
                            if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                dtf = -1
                                break
                        if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                            if Hr not in failed_interpolation.keys():
                                failed_interpolation[Hr] = ''
                            failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                            continue
                        else:
                            upper_limits[dt][Hr][farm] = 0.5*(float(upper_limits[dti][Hr][farm]) + float(upper_limits[dtf][Hr][farm]))
                elif Hr in upper_limits[dt].keys():
                    for farm in list_of_farms:
                        if farm not in upper_limits[dt][Hr].keys() or upper_limits[dt][Hr][farm] == '':
                            if Hr not in missing_upper_limits.keys():
                                missing_upper_limits[Hr]= ''
                            missing_upper_limits[Hr] = missing_upper_limits[Hr] + farm + ', '
                            dti, i = dt, 1
                            while True: # We find the closest lower neighbor
                                dti = dt - timedelta(minutes=i*60)
                                if dti in upper_limits.keys() and Hr in upper_limits[dti].keys() and farm in upper_limits[dti][Hr].keys() and upper_limits[dti][Hr][farm]!='':
                                    break
                                else:
                                    i = i+1
                                if i > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dti = -1
                                    break
                            dtf, f = dt, 1
                            while True: # We find the closest upper neighbor
                                dtf = dt + timedelta(minutes=f*60)
                                if dtf in upper_limits.keys() and Hr in upper_limits[dtf].keys() and farm in upper_limits[dtf][Hr].keys() and upper_limits[dtf][Hr][farm]!='':
                                    break
                                else:
                                    f = f+1
                                if f > options.max_nb_interpolated_hours: # if we don't have data for 24 hours there is a bigger issue
                                    dtf = -1
                                    break
                            if dti == -1 or dtf == -1 or (i+f-1 > options.max_nb_interpolated_hours):
                                if Hr not in failed_interpolation.keys():
                                    failed_interpolation[Hr] = ''
                                failed_interpolation[Hr] = failed_interpolation[Hr] + farm + ', '
                                continue
                            else:
                                upper_limits[dt][Hr][farm] = 0.5*(float(upper_limits[dti][Hr][farm]) + float(upper_limits[dtf][Hr][farm]))
            visited_hr = []
            for Hr1 in missing_upper_limits.keys():
                if Hr1 in visited_hr:
                    continue
                horizons = [Hr1]
                for Hr2 in missing_upper_limits.keys():
                    if Hr2 != Hr1 and missing_upper_limits[Hr2] == missing_upper_limits[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print('***WARNING: Upper limit is missing for date ' + str(dt.date()) +  ' and time (UTC) ' + str(dt.time()) + ' for forecast horizon(s) ' +  hor.strip(', ') +  ' for farm(s) ' + missing_upper_limits[Hr1].strip(', ') + '. Interpolating using one hour increments...')
                visited_hr = visited_hr + horizons
            visited_hr_fi = []
            for Hr1 in failed_interpolation.keys():
                if Hr1 in visited_hr_fi:
                    continue
                horizons = [Hr1]
                for Hr2 in failed_interpolation.keys():
                    if Hr2 != Hr1 and failed_interpolation[Hr2] == failed_interpolation[Hr1]:
                        horizons.append(Hr2)
                horizons = sorted(horizons)
                if len(horizons) == 1:
                    hor = str(horizons[0])
                else:
                    h0 = horizons[0]
                    h1 = horizons[0]
                    hor = ''
                    for i in range(1, len(horizons)):
                        if horizons[i] == h1+1:
                            h1 = horizons[i]
                            if i != len(horizons)-1:
                                continue
                        if h0 == h1:
                            hor = hor + str(h0) + ', '
                        elif h1 == h0+1:
                            hor = hor + str(h0) + ', ' + str(h1) + ', '
                        else:
                            hor = hor + str(h0) + ' - ' + str(h1) + ', '
                        h0 = horizons[i]
                        h1 = horizons[i]
                print("***WARNING: Interpolation failed - we can't find neither lower nor upper neighbors for date time (UTC) " + str(dt) + ", within " + str(options.max_nb_interpolated_hours) + " hours, for forecast horizon(s) " + hor.strip(', ') + " and farms: " + failed_interpolation[Hr1].strip(', ') + ".")
                visited_hr_fi = visited_hr_fi + horizons
                    
print('')                
                                             
#####################
##### aggregate #####
#####################
print('Aggregating..')
aggregated_actuals = OrderedDict()  # {date: actuals for aggregated farms}
for dt in actuals.keys():
    actuals_aggregated_farms = 0
    flag = 1
    missing_actuals = ''
    for farm in list_of_farms:
        if farm not in actuals[dt].keys() or actuals[dt][farm] == '' or actuals[dt][farm] == 'Shutdown':
            missing_actuals = missing_actuals + farm + ', ' 
            flag = 0
        if flag == 1:
            actuals_aggregated_farms = actuals_aggregated_farms + float(actuals[dt][farm])
    if flag == 0:
        aggregated_actuals[dt] = ''
        print('Deleting date time', dt, 'from actuals due to missing actuals for farms:', missing_actuals.strip(', '))
    else:
        aggregated_actuals[dt] = actuals_aggregated_farms

aggregated_means = OrderedDict()    # {date when forecast was made: {forecast horizon: mean for aggregated farms}}
for dt in forecast.keys():
    missing_forecasts = OrderedDict()
    if Hr not in forecast[dt].keys():
        missing_forecasts[Hr] = 'all farms'
    else:
        for Hr in forecast[dt].keys():
            mean_aggregated_farms = 0
            flag = 1
            for farm in list_of_farms:
                if farm not in forecast[dt][Hr].keys() or forecast[dt][Hr][farm] == '':
                    if Hr not in missing_forecasts.keys():
                        missing_forecasts[Hr] = ''
                    missing_forecasts[Hr] = missing_forecasts[Hr] + farm + ', '
                    flag = 0
                if flag ==1:
                    mean_aggregated_farms = mean_aggregated_farms + float(forecast[dt][Hr][farm])
            if dt not in aggregated_means.keys():
                aggregated_means[dt] = OrderedDict()
            if flag == 0:
                aggregated_means[dt][Hr] = ''
            else:
                aggregated_means[dt][Hr] = mean_aggregated_farms
    visited_hr = []
    for Hr1 in missing_forecasts.keys():
        if Hr1 in visited_hr:
            continue
        horizons = [Hr1]
        for Hr2 in missing_forecasts.keys():
            if Hr2 != Hr1 and missing_forecasts[Hr2] == missing_forecasts[Hr1]:
                horizons.append(Hr2)
        horizons = sorted(horizons)
        if len(horizons) == 1:
            hor = str(horizons[0])
        else:
            h0 = horizons[0]
            h1 = horizons[0]
            hor = ''
            for i in range(1, len(horizons)):
                if horizons[i]==h1+1:
                    h1 = horizons[i]
                    if i != len(horizons)-1:
                        continue
                if h0 == h1:
                    hor = hor + str(h0) + ', '
                elif h1 == h0+1:
                    hor = hor + str(h0) + ', ' + str(h1) + ', '
                else:
                    hor = hor + str(h0) + ' - ' + str(h1) + ', '
                h0 = horizons[i]
                h1 = horizons[i]
        print('Deleting date time', dt, 'and forecast horizon(s)', hor.strip(', '), 'from forecasts due to missing forecasts for farms:', missing_forecasts[Hr1].strip(', '))
        visited_hr = visited_hr + horizons
        
aggregated_lower_limits = OrderedDict()    # {date when forecast was made: {forecast horizon: lower limit for aggregated farms}}
for dt in lower_limits.keys():
    missing_lower_limits = OrderedDict()
    if Hr not in lower_limits[dt].keys():
        missing_lower_limits[Hr] = 'all farms'
    else:
        for Hr in lower_limits[dt].keys():
            lower_limits_aggregated_farms = 0
            flag = 1
            for farm in list_of_farms:
                if farm not in lower_limits[dt][Hr].keys() or lower_limits[dt][Hr][farm] == '':
                    if Hr not in missing_lower_limits.keys():
                        missing_lower_limits[Hr] = ''
                    missing_lower_limits[Hr] = missing_lower_limits[Hr] + farm + ', '
                    flag = 0
                if flag == 1:
                    lower_limits_aggregated_farms = lower_limits_aggregated_farms + float(lower_limits[dt][Hr][farm])
            if dt not in aggregated_lower_limits.keys():
                aggregated_lower_limits[dt] = OrderedDict()
            if flag == 0:
                aggregated_lower_limits[dt][Hr] = ''
            else:
                aggregated_lower_limits[dt][Hr] = lower_limits_aggregated_farms
    visited_hr = []
    for Hr1 in missing_lower_limits.keys():
        if Hr1 in visited_hr:
            continue
        horizons = [Hr1]
        for Hr2 in missing_lower_limits.keys():
            if Hr2 != Hr1 and missing_lower_limits[Hr2] == missing_lower_limits[Hr1]:
                horizons.append(Hr2)
        horizons = sorted(horizons)
        if len(horizons) == 1:
            hor = str(horizons[0])
        else:
            h0 = horizons[0]
            h1 = h0
            hor = ''
            for i in range(1, len(horizons)):
                if horizons[i]==h1+1:
                    h1 = horizons[i]
                    if i != len(horizons)-1:
                        continue
                if h0 == h1:
                    hor = hor + str(h0) + ', '
                elif h1 == h0+1:
                    hor = hor + str(h0) + ', ' + str(h1) + ', '
                else:
                    hor = hor + str(h0) + ' - ' + str(h1) + ', '
                h0 = horizons[i]
                h1 = horizons[i]
        print('Deleting date time', dt, 'and forecast horizon(s)', hor.strip(', '), 'from lower limits due to missing lower limits for farms:', missing_lower_limits[Hr1].strip(', '))
        visited_hr = visited_hr + horizons

aggregated_upper_limits = OrderedDict()    # {date when forecast was made: {forecast horizon: upper limit for aggregated farms}}
for dt in upper_limits.keys():
    missing_upper_limits = OrderedDict()
    if Hr not in upper_limits[dt].keys():
        missing_upper_limits[Hr] = 'all farms'
    else:
        for Hr in upper_limits[dt].keys():
            upper_limits_aggregated_farms = 0
            flag = 1
            for farm in list_of_farms:
                if farm not in upper_limits[dt][Hr].keys() or upper_limits[dt][Hr][farm] == '':
                    if Hr not in missing_upper_limits.keys():
                        missing_upper_limits[Hr] = ''
                    missing_upper_limits[Hr] = missing_upper_limits[Hr] + farm + ', '
                    flag = 0
                if flag == 1:
                    upper_limits_aggregated_farms = upper_limits_aggregated_farms + float(upper_limits[dt][Hr][farm])
            if dt not in aggregated_upper_limits.keys():
                aggregated_upper_limits[dt] = OrderedDict()
            if flag == 0:
                aggregated_upper_limits[dt][Hr] = ''
            else:
                aggregated_upper_limits[dt][Hr] = upper_limits_aggregated_farms
    visited_hr = []
    for Hr1 in missing_upper_limits.keys():
        if Hr1 in visited_hr:
            continue
        horizons = [Hr1]
        for Hr2 in missing_upper_limits.keys():
            if Hr2 != Hr1 and missing_upper_limits[Hr2] == missing_upper_limits[Hr1]:
                horizons.append(Hr2)
        horizons = sorted(horizons)
        if len(horizons) == 1:
            hor = str(horizons[0])
        else:
            h0 = horizons[0]
            h1 = h0
            hor = ''
            for i in range(1, len(horizons)):
                if horizons[i]==h1+1:
                    h1 = horizons[i]
                    if i != len(horizons)-1:
                        continue
                if h0 == h1:
                    hor = hor + str(h0) + ', '
                elif h1 == h0+1:
                    hor = hor + str(h0) + ', ' + str(h1) + ', '
                else:
                    hor = hor + str(h0) + ' - ' + str(h1) + ', '
                h0 = horizons[i]
                h1 = horizons[i]
        print('Deleting date time', dt, 'and forecast horizon(s)', hor.strip(', '), 'from upper limits due to missing upper limits for farms:', missing_upper_limits[Hr1].strip(', '))
        visited_hr = visited_hr + horizons
        

#########################
##### write to file #####
#########################
if options.actuals_writing_file:
    print('Writing actuals to file..')
    if not os.path.exists(options.actuals_output_directory): 
        os.mkdir(options.actuals_output_directory)
    output_directory = options.actuals_output_directory + os.sep
    file1 = open(output_directory + options.actuals_writing_file, 'w')
    file1.write('#dt,,observed value' + '\n')
    for dt in aggregated_actuals.keys():
        date = dt.strftime('%m/%d/%y %H:%M')
        file1.write(date + ',,' + str(aggregated_actuals[dt]) + '\n')
    file1.close()
    

if options.forecast_writing_file:
    print('Writing forecasts to file..')
    if not os.path.exists(options.forecast_output_directory): 
        os.mkdir(options.forecast_output_directory)
    output_directory = options.forecast_output_directory + os.sep
    file2 = open(output_directory + options.forecast_writing_file, 'w')
    file2.write('# BPA WIT Initiative: Centralized Wind Power Forecasting,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Aggregated wind power forecasts for all wind generation in the BPA BAA,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Historical data set,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Data Description: estimated --> AVERAGE <-- generation,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Notes:,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Date/Time is in UTC or Greenwich Mean Time and is the time the forecast was produced,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# HR01 is the first hour of the forecast or next hour,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# Inquiries at BPA Wind Integration Team / BPAwindintegration@bpa.gov,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('# ----------------------------------------------------------------------------,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
    file2.write('Date/Time (UTC),Hr01,Hr02,Hr03,Hr04,Hr05,Hr06,Hr07,Hr08,Hr09,Hr10,Hr11,Hr12,Hr13,Hr14,Hr15,Hr16,Hr17,Hr18,Hr19,Hr20,Hr21,Hr22,Hr23,Hr24,Hr25,Hr26,Hr27,Hr28,Hr29,Hr30,Hr31,Hr32,Hr33,Hr34,Hr35,Hr36,Hr37,Hr38,Hr39,Hr40,Hr41,Hr42,Hr43,Hr44,Hr45,Hr46,Hr47,Hr48,Hr49,Hr50,Hr51,Hr52,Hr53,Hr54,Hr55,Hr56,Hr57,Hr58,Hr59,Hr60,Hr61,Hr62,Hr63,Hr64,Hr65,Hr66,Hr67,Hr68,Hr69,Hr70,Hr71,Hr72\n')
    for dt in aggregated_means.keys():
        date = dt.strftime('%m/%d/%Y %H:%M')
        string_of_means = ''
        for Hr in aggregated_means[dt].keys():
            string_of_means = string_of_means + ',' + str(aggregated_means[dt][Hr]).split('.')[0]
        file2.write(date + string_of_means + '\n')
    file2.close()
    
if options.lower_limit_writing_file:
    print('Writing lower limtis to file..')
    if not os.path.exists(options.lower_limit_output_directory): 
        os.mkdir(options.lower_limit_output_directory)
    output_directory = options.lower_limit_output_directory + os.sep
    file3 = open(output_directory + options.lower_limit_writing_file, 'w')
    file3.write('Date/Time (UTC),Hr01,Hr02,Hr03,Hr04,Hr05,Hr06,Hr07,Hr08,Hr09,Hr10,Hr11,Hr12,Hr13,Hr14,Hr15,Hr16,Hr17,Hr18,Hr19,Hr20,Hr21,Hr22,Hr23,Hr24,Hr25,Hr26,Hr27,Hr28,Hr29,Hr30,Hr31,Hr32,Hr33,Hr34,Hr35,Hr36,Hr37,Hr38,Hr39,Hr40,Hr41,Hr42,Hr43,Hr44,Hr45,Hr46,Hr47,Hr48,Hr49,Hr50,Hr51,Hr52,Hr53,Hr54,Hr55,Hr56,Hr57,Hr58,Hr59,Hr60,Hr61,Hr62,Hr63,Hr64,Hr65,Hr66,Hr67,Hr68,Hr69,Hr70,Hr71,Hr72\n')
    for dt in aggregated_lower_limits.keys():
        date = dt.strftime('%m/%d/%Y %H:%M')
        string_of_lower_limits = ''
        for Hr in aggregated_lower_limits[dt].keys():
            string_of_lower_limits = string_of_lower_limits + ',' + str(aggregated_lower_limits[dt][Hr]).split('.')[0]
        file3.write(date + string_of_lower_limits + '\n')
    file3.close()
    
if options.upper_limit_writing_file:
    print('Writing upper limits to file..')
    if not os.path.exists(options.upper_limit_output_directory): 
        os.mkdir(options.upper_limit_output_directory)
    output_directory = options.upper_limit_output_directory + os.sep
    file4 = open(output_directory + options.upper_limit_writing_file, 'w')
    file4.write('Date/Time (UTC),Hr01,Hr02,Hr03,Hr04,Hr05,Hr06,Hr07,Hr08,Hr09,Hr10,Hr11,Hr12,Hr13,Hr14,Hr15,Hr16,Hr17,Hr18,Hr19,Hr20,Hr21,Hr22,Hr23,Hr24,Hr25,Hr26,Hr27,Hr28,Hr29,Hr30,Hr31,Hr32,Hr33,Hr34,Hr35,Hr36,Hr37,Hr38,Hr39,Hr40,Hr41,Hr42,Hr43,Hr44,Hr45,Hr46,Hr47,Hr48,Hr49,Hr50,Hr51,Hr52,Hr53,Hr54,Hr55,Hr56,Hr57,Hr58,Hr59,Hr60,Hr61,Hr62,Hr63,Hr64,Hr65,Hr66,Hr67,Hr68,Hr69,Hr70,Hr71,Hr72\n')
    for dt in aggregated_upper_limits.keys():
        date = dt.strftime('%m/%d/%Y %H:%M')
        string_of_upper_limits = ''
        for Hr in aggregated_upper_limits[dt].keys():
            string_of_upper_limits = string_of_upper_limits + ',' + str(aggregated_upper_limits[dt][Hr]).split('.')[0]
        file4.write(date + string_of_upper_limits + '\n')
    file4.close()
            
        
