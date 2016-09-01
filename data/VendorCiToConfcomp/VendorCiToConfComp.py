

import sys, os, datetime, time, numpy
import VendorCiToConfcompOptions

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
    input_file_list_lower_limits = []
    input_file_list_upper_limits = []
    
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
        
        # read the type of the file - ll = lower limit, ul = upper limit
        type_of_file = line.split(' ')[0]
        
        # check if file already in list of input files
        if name_of_file in input_file_list_lower_limits or name_of_file in input_file_list_upper_limits:
            continue
        
        # add files and their directory to the according list
        if type_of_file == 'll':
            input_file_list_lower_limits.append(name_of_file)
        elif type_of_file == 'ul':
            input_file_list_upper_limits.append(name_of_file)
    
    # close the file containing the list of input files
    list_file.close()

    # return the list of input files
    return input_file_list_lower_limits, input_file_list_upper_limits


###
# Read the Chosen input files and store the data in a dictionary
###
def get_unprocessed_ci_limits(options, list_of_input_files_lower_limits, list_of_input_files_upper_limits):
    print('read unprocessed lower limits..')
    
    # create empty dictionary where the unprocessed lower limit will be stored
    unprocessed_lower_limit = OrderedDict()

    # Test if the user has given a reasonable input for the option of where the lower limit has to be taken.
    # Change the format from int to datetime
    try:
        time_of_lower_limit = datetime.datetime.strptime(str(options.time_of_forecast), '%H')
    # raise error if the given input hasn't been an integer between 0 and 23
    except:
        print("***ERROR: --time-of-forecast has to be a number in digits between 0 and 23")
        sys.exit(1)
    
    
    # read the files in the list of input files 
    for file_name in list_of_input_files_lower_limits:
        
        # get the directory of the input files
        input_file_directory = options.input_directory + os.sep # by default: 'input_files'
        
        # open the files in the list of input files
        f = open(input_file_directory + file_name,'r')
        
        # read the files in the list of input files
        for lower_limit_line in f:
            
            # ignore lines of no interest
            if not lower_limit_line[0].isdigit():
                continue
                
            # store the date + time
            dt = datetime.datetime.strptime(lower_limit_line.split(',')[0], '%m/%d/%Y %H:%M')
            
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
            if not dt.hour == time_of_lower_limit.hour:
                continue
            
            print('dt UTC for lower limit =', dt+datetime.timedelta(hours=difference_pdt_to_utc))
            print('dt PDT for lower limit=', dt)
            print(lower_limit_line)
            
            # search in the row of lower limit values for the beginning of the next day.
            lower_limit_place_of_hour_0 = 24 - time_of_lower_limit.hour
            
            # populate dictionary with lower limits
            if options.intraday_forecast :
                 for i in range(1,23-options.time_of_forecast+1):
                    # check if lower limit is in file and react if data is missing
                    if lower_limit_line.split(',')[i] == '' or not lower_limit_line.split(',')[i][0].isdigit():
                        continue
                    
                    # save the lower limit of the next HOURS after time forecast of this day, this is 1 based hours (or periods in the day)
                    next_hour= dt+datetime.timedelta(hours=i)
                    print('date predicted =', next_hour)
                    unprocessed_lower_limit[next_hour] = float(lower_limit_line.split(',')[i])
                    print('lower limit = ', unprocessed_lower_limit[next_hour])
            else :
                for i in range(24):
                    # check if lower limit is in file and react if data is missing
                    if lower_limit_line.split(',')[lower_limit_place_of_hour_0+i] == '' or not lower_limit_line.split(',')[lower_limit_place_of_hour_0+i][0].isdigit():
                        continue
                    
                    # save the lower limit of the next day
                    hour_of_next_day = dt+datetime.timedelta(hours=lower_limit_place_of_hour_0+i)
                    print('date predicted =', hour_of_next_day)
                    unprocessed_lower_limit[hour_of_next_day] = float(lower_limit_line.split(',')[lower_limit_place_of_hour_0+i])
                    print('lower limit = ', unprocessed_lower_limit[hour_of_next_day])
                
        f.close()
        
    print('read unprocessed upper limits..')
    
    # create empty dictionary where the unprocessed upper limit will be stored
    unprocessed_upper_limit = OrderedDict()

    # Test if the user has given a reasonable input for the option of where the upper limit has to be taken.
    # Change the format from int to datetime
    try:
        time_of_upper_limit = datetime.datetime.strptime(str(options.time_of_forecast), '%H')
    # raise error if the given input hasn't been an integer between 0 and 23
    except:
        print("***ERROR: --time-of-forecast has to be a number in digits between 0 and 23")
        sys.exit(1)
    
    
    # read the files in the list of input files 
    for file_name in list_of_input_files_upper_limits:
        
        # get the directory of the input files
        input_file_directory = options.input_directory + os.sep # by default: 'input_files'
        
        # open the files in the list of input files
        f = open(input_file_directory + file_name,'r')
        
        # read the files in the list of input files
        for upper_limit_line in f:
            
            # ignore lines of no interest
            if not upper_limit_line[0].isdigit():
                continue
                
            # store the date + time
            dt = datetime.datetime.strptime(upper_limit_line.split(',')[0], '%m/%d/%Y %H:%M')
            
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
            if not dt.hour == time_of_upper_limit.hour:
                continue
            
            print('dt UTC for upper limit =', dt+datetime.timedelta(hours=difference_pdt_to_utc))
            print('dt PDT for upper limit=', dt)
            print(upper_limit_line)
            
            # search in the row of upper limit values for the beginning of the next day.
            upper_limit_place_of_hour_0 = 24 - time_of_upper_limit.hour
            
            # populate dictionary with upper limits
            if options.intraday_forecast :
                 for i in range(1,23-options.time_of_forecast+1):
                    # check if upper limit is in file and react if data is missing
                    if upper_limit_line.split(',')[i] == '' or not upper_limit_line.split(',')[i][0].isdigit():
                        continue
                    
                    # save the upper limit of the next HOURS after time forecast of this day, this is 1 based hours (or periods in the day)
                    next_hour= dt+datetime.timedelta(hours=i)
                    print('date predicted =', next_hour)
                    unprocessed_upper_limit[next_hour] = float(upper_limit_line.split(',')[i])
                    print('upper limit = ', unprocessed_upper_limit[next_hour])
            else :
                for i in range(24):
                    # check if upper limit is in file and react if data is missing
                    if upper_limit_line.split(',')[upper_limit_place_of_hour_0+i] == '' or not upper_limit_line.split(',')[upper_limit_place_of_hour_0+i][0].isdigit():
                        continue
                    
                    # save the upper limit of the next day
                    hour_of_next_day = dt+datetime.timedelta(hours=upper_limit_place_of_hour_0+i)
                    print('date predicted =', hour_of_next_day)
                    unprocessed_upper_limit[hour_of_next_day] = float(upper_limit_line.split(',')[upper_limit_place_of_hour_0+i])
                    print('upper limit = ', unprocessed_upper_limit[hour_of_next_day])
                
        f.close()
    return unprocessed_lower_limit, unprocessed_upper_limit


###
# Write the data into a file
###
def write_limit_file(options, lower_limit_dictionary, upper_limit_dictionary):
    print('write lower limits into file..')
    
    print(options.lower_limit_output_directory)
    # create directory for output files given by user if it does not exist
    if not os.path.exists(options.lower_limit_output_directory): 
        os.mkdir(options.lower_limit_output_directory)
    
    # open/create file to write into
    output_directory = options.lower_limit_output_directory + os.sep
    name_of_lower_limit_file = options.name_of_lower_limit_file
    outputfile = open(output_directory + name_of_lower_limit_file,'w')
    
    # write description at the beginning of the file
    print('# Lower limit file written by ForecastToConfComp.py', file=outputfile)
    print('# time , lower limit', file=outputfile)
    print('#', file=outputfile)
    
    # write the lower limits into the file
    for dt in lower_limit_dictionary:
        
        # save the parts of the datetime as strings to stay with the output structure closer to the structure of the BPA files
        year = str(dt.year)
        year = year[2:]
        month = str(dt).split('-')[1]
        day = str(dt).split('-')[2].split(' ')[0]
        hour = str(dt).split(' ')[1].split(':')[0]
        minute = str(dt).split(' ')[1].split(':')[1]
        
        # save lower limits at dt as string
        lower_limit = str(lower_limit_dictionary[dt])
        print(month+'/'+day+'/'+year+' '+hour+':'+minute+','+lower_limit, file=outputfile)
        
    print('write upper limits into file..')
    
    print(options.upper_limit_output_directory)
    # create directory for output files given by user if it does not exist
    if not os.path.exists(options.upper_limit_output_directory): 
        os.mkdir(options.upper_limit_output_directory)
    
    # open/create file to write into
    output_directory = options.upper_limit_output_directory + os.sep
    name_of_upper_limit_file = options.name_of_upper_limit_file
    outputfile = open(output_directory + name_of_upper_limit_file,'w')
    
    # write description at the beginning of the file
    print('# Upper limit file written by VendorCiToConfComp.py', file=outputfile)
    print('# time , upper limit', file=outputfile)
    print('#', file=outputfile)
    
    # write the upper limits into the file
    for dt in upper_limit_dictionary:
        
        # save the parts of the datetime as strings to stay with the output structure closer to the structure of the BPA files
        year = str(dt.year)
        year = year[2:]
        month = str(dt).split('-')[1]
        day = str(dt).split('-')[2].split(' ')[0]
        hour = str(dt).split(' ')[1].split(':')[0]
        minute = str(dt).split(' ')[1].split(':')[1]
        
        # save lower limits at dt as string
        upper_limit = str(upper_limit_dictionary[dt])
        print(month+'/'+day+'/'+year+' '+hour+':'+minute+','+upper_limit, file=outputfile)



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
        options_parser = VendorCiToConfcompOptions.construct_options_parser()
        (options, args) = options_parser.parse_args(args=sys.argv)
    except SystemExit:
        # the parser throws a system exit if "-h" is specified - catch
        # it to exit gracefully.
        sys.exit(0)
    
    
    ###
    # Test necessary user input
    ###
    if options.name_of_lower_limit_file is None:
        print("***ERROR: The name of the lower limit output file has to be given.")
        print("          Specify the name of the lower limit output file with the option --name-of-lower-limit-file")
        sys.exit(1)
    if options.name_of_upper_limit_file is None:
        print("***ERROR: The name of the upper limit output file has to be given.")
        print("          Specify the name of the upper limit output file with the option --name-of-upper-limit-file")
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
ListOfInputFilesLowerLimits, ListOfInputFilesUpperLimits = get_list_of_files(options)
print(ListOfInputFilesUpperLimits)
###
# Reading the chosen input files and store the data in the dictionaries UnprocessedLowerLimit and UnprocessedUpperLimit
# it gets the limits for the 24 hours of the next day, made on the hour specified in the option the day before
###
UnprocessedLowerLimit, UnprocessedUpperLimit = get_unprocessed_ci_limits(options,ListOfInputFilesLowerLimits, ListOfInputFilesUpperLimits)

###
# writing the data into a new file
###
write_limit_file(options, UnprocessedLowerLimit, UnprocessedUpperLimit)

print('Done!')
