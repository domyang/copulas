

import sys, os, datetime, time, numpy
import ForecastToConfcompOptions

try:
    from collections import OrderedDict
except:
    from ordereddict import OrderedDict

# This program reads the given Forecast files and escludes days according to parameters given by the user.

########################################################################################################################
#                                               auxiliary functions                                                    #
########################################################################################################################
###
# Get daylight saving time (the datetime of that year)
###
# rule: 2nd Sunday in March
def get_daylight_saving_time(year):
    dls_day = datetime.datetime.strptime('3/8/'+year, '%m/%d/%Y')
    for i in range(20):
        if dls_day.month == 3 and dls_day.weekday() == 6 and week_of_month(dls_day) == 2:
            break
        dls_day += datetime.timedelta(days=1)
    dls_day += datetime.timedelta(hours=9) 
    return dls_day
    
###
# Get end of daylight saving (the datetime of that year)
###
# rule: 1st Sunday in November
def get_daylight_saving_end(year):
    dls_end = datetime.datetime.strptime('11/1/'+year, '%m/%d/%Y')
    for i in range(20):
        if dls_end.month == 11 and dls_end.weekday() == 6 and week_of_month(dls_end) == 1:
            break
        dls_end += datetime.timedelta(days=1)
    dls_end += datetime.timedelta(hours=8) 
    return dls_end

###
# Get the week of the month of the given date
###
def week_of_month(date):
    month = date.month
    week = 0
    while date.month == month:
        week += 1
        date -= datetime.timedelta(days=7)
    return week



########################################################################################################################
#                                                  major functions                                                     #
########################################################################################################################
###
# Read the list of files to be read
###
def get_list_of_files(options):
    input_file_list = []
    process_file_list = []
    # get the file of the list with the input files
    list_file_directory = options.list_of_input_file_name #by deafault 'list_of_input_files.dat'
    
    # open the file containing the list of input files
    list_file = open(list_file_directory,'r')
    
    # read the file containing the list of input files
    for line in list_file:
        # ignore lines commented out
        if line[0] == '#':
            continue
        
        # read the name of the file
        name_of_file = str(line.split('\n')[0].split(' ')[1])
        
        # read the type of the file - fc = forecast, br = balancing reserves
        type_of_file = line.split(' ')[0]
        
        # check if file already in list of input files
        if name_of_file in input_file_list:
            continue
        
        # add files and their directory to the according list
        if type_of_file == 'fc':
            input_file_list.append(name_of_file)
        elif type_of_file == 'br':
            process_file_list.append(name_of_file)
    
    # close the file containing the list of input files
    list_file.close()

    # return the list of input files
    return input_file_list, process_file_list


###
# Read the Chosen input files and store the data in a dictionary
###
def get_unprocessed_forecast(options,list_of_input_files):
    print('read unprocessed forecast..')
    
    # create empty dictionary where the unprocessed forecast will be stored
    unprocessed_forecast = OrderedDict()

    # Test if the user has given a reasonable input for the option of where the forecast has to be taken.
    # Change the format from int to datetime
    try:
        time_of_forecast = datetime.datetime.strptime(str(options.time_of_forecast), '%H')
    # raise error if the given input hasn't been an integer between 0 and 23
    except:
        print("***ERROR: --time-of-forecast has to be a number in digits between 0 and 23")
        sys.exit(1)
    
    
    # read the files in the list of input files 
    for file_name in list_of_input_files:
        
        # get the directory of the input files
        input_file_directory = options.input_directory + os.sep # by default: 'input_files'
        
        # open the files in the list of input files
        f = open(input_file_directory + file_name,'r')
        
        # read the files in the list of input files
        for forecast_line in f:
            
            # ignore lines of no interest
            if not forecast_line[0].isdigit():
                continue
                
            # store the date + time
            dt = datetime.datetime.strptime(forecast_line.split(',')[0], '%m/%d/%Y %H:%M')
            
            # get the day and time of the start daylight saving in the current year
            daylight_saving_start = get_daylight_saving_time(str(dt.year))
            
            # get the day and time of the end daylight saving in the current year
            daylight_saving_end = get_daylight_saving_end(str(dt.year))


            # get the difference between pacific time and utc
            difference_pdt_to_utc = 'will be set now' # with respect to daylight saving
            if dt >= daylight_saving_start and dt <= daylight_saving_end:
                difference_pdt_to_utc = 7
            elif dt < daylight_saving_start or dt > daylight_saving_end:
                difference_pdt_to_utc = 8
            
            # CHANGES THE TIME FROM GMT TO PDT! NOT TO THE HOUR OF THE FORECAST!
            dt -= datetime.timedelta(hours=difference_pdt_to_utc)
            
            # search in the document for time where the forecast should be made. 
            # i.e.: skip every line where the forecast was made at a different time.
            if not dt.hour == time_of_forecast.hour:
                continue
            
            print('dt UTC for forecast =', dt+datetime.timedelta(hours=difference_pdt_to_utc))
            print('dt PDT for forecast=', dt)
            print(forecast_line)
            
            # search in the row of forecast values for the beginning of the next day.
            forecast_place_of_hour_0 = 24 - time_of_forecast.hour
            
            # populate dictionary with forecast
            if options.intraday_forecast :
                 for i in range(1,23-options.time_of_forecast+1):
                    # check if forecast is in file and react if data is missing
                    if forecast_line.split(',')[i] == '' or  not forecast_line.split(',')[i][0].isdigit():
                        continue
                    
                    # save the forecast of the next HOURS after time forecast of this day, this is 1 based hours (or periods in the day)
                    next_hour= dt+datetime.timedelta(hours=i)
                    print('date predicted =', next_hour)
                    unprocessed_forecast[next_hour] = float(forecast_line.split(',')[i])
                    print('forecast = ', unprocessed_forecast[next_hour])
            else :
                for i in range(24):
                    # check if forecast is in file and react if data is missing
                    if forecast_line.split(',')[forecast_place_of_hour_0+i] == '' or not forecast_line.split(',')[forecast_place_of_hour_0+i][0].isdigit():
                        continue
                    
                    # save the forecast of the next day
                    hour_of_next_day = dt+datetime.timedelta(hours=forecast_place_of_hour_0+i)
                    print('date predicted =', hour_of_next_day)
                    unprocessed_forecast[hour_of_next_day] = float(forecast_line.split(',')[forecast_place_of_hour_0+i])
                    print('forecast = ', unprocessed_forecast[hour_of_next_day])
                
        f.close()
    return unprocessed_forecast


###
# Process the forecast
###
def process_forecast(options, unprocessed_forecast, process_file_list):
    print('process the forecast..')
    
    # generate necessary dictionaries
    anomalistic_times = {}
    days_with_anomalies = {}
    
    # read the options to process the forecast
    critical_number_of_times = options.critical_number_of_times
    negative_boundary        = float(options.negative_boundary)
    positive_boundary        = options.positive_boundary
    
    # directory for process data
    file_directory = 'balancing_reserves' + os.sep
    # open data with information to process the forecast
    for file_name in process_file_list:
    
        # open the file with the information
        f = open(file_directory + file_name,'r')
    
        # read the files in the list of the process files
        for process_line in f:
        
            # ignore lines of no interest
            if not process_line[0].isdigit():
                continue
            # ignore the table that explains the numbers in the file
            if len(process_line.split('/')) == 1:
                continue
            
            # ignore dates where data is missing
            #wind_state
            if process_line.split(',')[4] == '':
                continue
                
            # save the data into variables
            date                  = datetime.datetime.strptime(process_line.split(',')[0], '%m/%d/%y %H:%M')
            wind_state            = float(process_line.split(',')[4])
            if wind_state <= negative_boundary or wind_state >= positive_boundary:
                anomalistic_times[date] = wind_state

    # create a dictionary, whith the datetime at 00:00 are the keys
    # For every day the value is a list of all times with anomalies
    for dt_ano in anomalistic_times:
        # save the day the anomalistic time is in
        day_of_time = dt_ano
        # set it to 00:00 so we won't test in several times
        day_of_time -= datetime.timedelta(hours=dt_ano.hour)
        day_of_time -= datetime.timedelta(minutes=dt_ano.minute) 

        if day_of_time not in days_with_anomalies:
            days_with_anomalies[day_of_time] = [dt_ano]
        else:
            days_with_anomalies[day_of_time].append(dt_ano)
    
    # sort out all days which need to be sorted out
    # The number of anomalies in a day necessary until we exclude a day is given with --
    for dt_at0 in list(days_with_anomalies.keys()):
        if len(days_with_anomalies[dt_at0]) < critical_number_of_times:
            del days_with_anomalies[dt_at0]

    # delete the days with too many anomalities out of the forecast dictionary
    for dt_fc in list(unprocessed_forecast.keys()):
        for dt_at0 in days_with_anomalies:
            if dt_fc.year == dt_at0.year and dt_fc.month == dt_at0.month and dt_fc.day == dt_at0.day:
                del unprocessed_forecast[dt_fc]

    return unprocessed_forecast
    

###
# Write the data into a file
###
def write_forecast_file(options, forecast_dictionary):
    print('write forecast into file..')
    
    print(options.output_directory)
    # create directory for output files given by user if it does not exist
    if not os.path.exists(options.output_directory): 
        os.mkdir(options.output_directory)
    
    # open/create file to write into
    output_directory = options.output_directory + os.sep
    name_of_forecast_file = options.name_of_forecast_file
    outputfile = open(output_directory + name_of_forecast_file,'w')
    
    # write description at the beginning of the file
    print('# Forecast file written by ForecastToConfComp.py', file=outputfile)
    print('# time , forecast', file=outputfile)
    print('#', file=outputfile)
    
    # write the forecast into the file
    for dt in forecast_dictionary:
        # save the parts of the datetime as strings to stay with the output structure closer to the structure of the BPA files
        year = str(dt.year)
        year = year[2:]
        month = str(dt).split('-')[1]
        day = str(dt).split('-')[2].split(' ')[0]
        hour = str(dt).split(' ')[1].split(':')[0]
        minute = str(dt).split(' ')[1].split(':')[1]
        
        # save forecast at dt as string
        forecast = str(forecast_dictionary[dt])
        print(month+'/'+day+'/'+year+' '+hour+':'+minute+','+forecast, file=outputfile)



########################################################################################################################
#                                                   main function                                                      #
########################################################################################################################
if __name__ == '__main__':

    ##############################
    # MAIN EXECUTION STARTS HERE #
    ##############################
    
    ###
    # parse the command-line options 
    ###
    try:
        options_parser = ForecastToConfcompOptions.construct_options_parser()
        (options, args) = options_parser.parse_args(args=sys.argv)
    except SystemExit:
        # the parser throws a system exit if "-h" is specified - catch
        # it to exit gracefully.
        sys.exit(0)
    
    
    ###
    # Test necessary user input
    ###
    if options.name_of_forecast_file is None:
        print("***ERROR: The name of the output file has to be given.")
        print("          Specify the name of the output file with the option --name-of-forecast-file")
        sys.exit(1)
    
###
# Get a list with all the files which should be processed by this program. 
###
# 
# These files have to be written in a list in the file specified by the option --list-of-input-file-name (which is by default 'list_of_input_files.dat') in the same directory as this program
# and stored in a folder specified by the option --input-directory 
# by default this folder has the name 'input_files'
#
###
ListOfInputFiles,ProcessFileList = get_list_of_files(options)

###
# Reading the chosen input files and store the data in the dictionary UnprocessedForecast
# it gets the forecast for the 24 hours of the next day, made on the hour specified in the option the day before
###
UnprocessedForecast = get_unprocessed_forecast(options,ListOfInputFiles)

###
# Process the data
# remove the days that have too many anomalies, that is, don t fit in the boundaries given in option
###
ProcessedForecast = process_forecast(options, UnprocessedForecast, ProcessFileList)

###
# writing the data into a new file
###
#
# for now I write the unprocessed data into the file specified by the user
#
###
write_forecast_file(options, UnprocessedForecast)



print('Done!')





