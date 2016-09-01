import time
import math
import sys
import os
import string
import shutil

sys.dont_write_bytecode = True

from optparse import OptionParser, OptionGroup

def construct_options_parser():
    parser = OptionParser()
    
    directory_options = OptionGroup(parser, "Directory Options")
    forecast_options = OptionGroup(parser, "Forecast Options")
    output_options = OptionGroup(parser, "Output Options")
    
    parser.add_option_group(directory_options)
    parser.add_option_group(forecast_options)
    parser.add_option_group(output_options)
    
    directory_options.add_option("--list-of-input-file-name",
                            help="The name of the file containing the list with the input files",
                            action="store",
                            dest="list_of_input_file_name",
                            type="string",
                            default="list_of_input_files.dat")
    
    directory_options.add_option("--input-directory",
                            help="The directory where the input data is stored.The default is: 'input_files'",
                            action="store",
                            dest="input_directory",
                            type="string",
                            default="input_files")

    directory_options.add_option("--lower-limit-output-directory",
                            help="The directory where the lower limit output data is stored.The default is: 'lower_limit_output_files'",
                            action="store",
                            dest="lower_limit_output_directory",
                            type="string",
                            default="lower_limit_output_files")
                
    directory_options.add_option("--upper-limit-output-directory",
                            help="The directory where the upper limit output data is stored.The default is: 'upper_limit_output_files'",
                            action="store",
                            dest="upper_limit_output_directory",
                            type="string",
                            default="upper_limit_output_files")
       
    forecast_options.add_option("--intraday-forecast",
                                 help="If specified, we will consider that the forecast is made the same day as the hours forecasted, so we will have only data for the hours ahead.",
                                 action="store_true",
                                 dest="intraday_forecast",
                                 default=False)

    forecast_options.add_option("--time-of-forecast",
                                 help="The time (as an hour of the day) the forecast should be computed. Has to be given as an int between 0 and 23. The Default is 0",
                                 action="store",
                                 dest="time_of_forecast",
                                 type="int",
                                 default=None)

    output_options.add_option("--name-of-lower-limit-file",
                                 help="The name of the lower limit output file.",
                                 action="store",
                                 dest="name_of_lower_limit_file",
                                 type="string",
                                 default=None)
                                 
    output_options.add_option("--name-of-upper-limit-file",
                                 help="The name of the upper limit output file.",
                                 action="store",
                                 dest="name_of_upper_limit_file",
                                 type="string",
                                 default=None)
                                

    return parser
    
if __name__ == '__main__':
    print("ForErrorOptions.py cannot run from the command line.")
    sys.exit(1)
