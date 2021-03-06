from datetime import datetime, timedelta

import dateutil.parser as dt
import matplotlib.pyplot as plt
import utilities as ut


# READ ME:
#
# In this class, we define a copula class specifically for ONE series of data and multiple offsets
# The goal is to take into account the dependencies between what happened at t0 and at t0+offset
# Say: wind power forecasts errors, at t=0 and t=5 hours
#
#
#

class CopulaManager:
    length = 0
    dim = 0
    vect = []  # a list of floats: the real data to be considered
    date = []  # a list of ints of ordered date
    data = {}  # a dictionary of lists, additional parameters considered during the creation of the window

    # dictionary defining 3 things:
    #                           the offsets (list)
    #                           the window (numerous arguments)
    #                           the kind of data that has been given in argument
    parameters = {"offsets": [], 'predicted_day': dt.parse('2000-01-01 00:00'),
                  'date_range': ('', ''), 'first_hour': (-11, 12), 'forecast': [], 'forecast_d': [],
                  'type': 'Wind', 'location': 'total', 'kind': 'error'}

    # dictionary stating the need for specific data for given parameters
    data_need = {'date_range': [], 'first_hour': [], 'forecast': ['forecast'], 'forecast_d': ['forecast_d']}

    # After we have resized the data with the window:
    unif = [[]]  # the uniformly distributed 'vect' -> copula points
    vectM = [[]]
    dateM = []
    dataM = {}
    lengthM = 0
    indexes = set()  # the indexes of 'vect' corresponding to the data for offset 0

    #######################
    ###                 ###
    ###     METHODS     ###
    ###                 ###
    #######################


    ### initialises the copula:
    # arguments:
    #   () date: a list (of int) of ordered dates
    #   () vect: a list of float: the real data to be considered
    #   () data: a dictionary of lists, additional parameters considered during the creation of the window
    #   () parameters: a dictionary defining 3 things:
    #                           the offsets (list)
    #                           the window (numerous arguments)
    #                           the kind of data that has been given in argument
    def __init__(self, dates, vect, data, parameters, var_function=None):

        if isinstance(dates[0], datetime):
            date = dates
        else:
            date = [dt.parse(date) for date in dates]

        #print('preparing the data')
        if len(date) == len(vect):

            vect_temp = []
            date_temp = []
            data_temp = {}
            for key in data:
                data_temp[key] = []

            vect_it = iter(vect)
            date_it = iter(date)
            data_it = {}
            for key in data:
                data_it[key] = iter(data[key])

            cur = date[0]
            while (cur < date[-1]):
                index = next(date_it)

                while cur < index:
                    vect_temp.append(None)
                    for key in data:
                        data_temp[key].append(None)
                    date_temp.append(cur)
                    cur += timedelta(hours=1)

                while cur > index:
                    print('Warning: possibly twice the same entry at date %s' % str(index))
                    index = next(date_it)
                    next(vect_it)
                    for key in data:
                        next(data_it[key])

                if cur == index:
                    date_temp.append(index)
                    vect_temp.append(next(vect_it))
                    for key in data:
                        data_temp[key].append(next(data_it[key]))

                cur += timedelta(hours=1)

            self.vect = vect_temp
            self.date = date_temp
            self.length = len(vect_temp)
            self.data = data_temp
            self.var_function = var_function

        else:
            raise (
            RuntimeError('vect and date should be same length, (here %d,%d respectively)' % (self.length, len(date))))

        self.update(parameters)

        # print('copula initialized')

    # the update functions update the window to match the parameters arguments
    # returns nothing
    # arguments:
    #   () parameters: a dictionary specifying the window
    #   () data: additional data to be taken into account
    def update(self, parameters):

        #print('checking the parameters arguments')
        data = self.data
        self.check_parameters(parameters, data)
        parameters = self.parameters
        self.define_window()

        # initializing the temporary variables
        dim = self.dim
        vectM = [[] for _ in range(dim)]
        indexes = []
        dataM = {}
        dateM = []

        for key in self.data:
            dataM[key] = [[] for _ in range(dim)]
        offsets = parameters['offsets']
        max_offset = max(offsets)

        # print('selecting the desired points')
        # selecting the desired points

        for i, date in enumerate(self.date):

            if i + max_offset >= self.length:
                break

            point = [self.vect[i + offset] for offset in offsets]  # The values of vectM
            # Dictionary of forecast and forecast_d
            more_data = {key: [data[key][i + offset] for offset in offsets] for key in self.data}


            if self.parameters['window'](point, date, more_data):

                for j, value in enumerate(point):
                    vectM[j].append(value)

                indexes.append(i)

                for key in more_data:
                    for j in range(dim):
                        dataM[key][j].append(more_data[key][j])
                dateM.append(date)

        self.dateM = dateM
        self.dataM = dataM
        self.lengthM = len(vectM[0])

        # if these are errors in Solar power forecasts, applying reverse variance transformation to get the correct values in vectM
        # (variances were scaled to 1 for each solar hour, now we multiply by the variance corresponding
        # to the solar hour of 'predicted_day' with offset)
        if (parameters['type'] == 'Solar') and (parameters['kind'] == 'error'):
            if ('hour_sol' in data) and (self.var_function is not None):
                if 'predicted_day' not in parameters:

                    if 'date_range' in parameters:
                        parameters['predicted_day'] = parameters['date_range'][0][1]
                    else:
                        parameters['predicted_day'] = dt.parse('2000-01-01 00:00')

                time = parameters['predicted_day']
                var_sol_temp = []
                rise, sete = ut.sun_rise_set('%d-%d-%d' % (time.year, time.month, time.day))

                for offset in offsets:
                    var_sol_temp.append(self.var_function(
                        ((offset + time.hour + time.minute / 60 + time.second / 3600) - rise) / (sete - rise) * 12))

                for i in range(dim):
                    vectM[i] = [x * var_sol_temp[i] for x in vectM[i]]

        self.vectM = vectM
        # creating the copula variables
        # print('creating the copula points')
        indexes = set(indexes)

        self.indexes = indexes
        self.unif = ut.uniforms(vectM, rand=False)

    ### pprint doesn't plot: it just list all the attributes:
    def pprint(self):

        print('\n### SCALAR PARAMETERS ### \n \n length= %d , dim= %d , lengthM= %d \n\n' % (
        self.length, self.dim, self.lengthM))
        print('### PRIMARY DATA ### \n \n-> vect (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
        len(self.vect), self.vect[0], self.vect[1], self.vect[2], self.vect[3], self.vect[4], self.vect[5]))
        print("-> date (len=%d) ['%s', '%s', '%s', ...]\n" % (
        len(self.date), self.date[0], self.date[1], self.date[2]))
        string = '-> data (len=%d) {' % len(self.data)
        for key in self.data.keys():
            string += "'%s', " % key
            string = string[:-2]
        string += '}\n'
        for key in self.data.keys():
            string += '   ' + key + ' (len=%d) [' % len(self.data[key])
            for i in range(min(6, len(self.data[key]))):
                string += '%0.2f, ' % self.data[key][i]
            string += '...]\n'
        print(string)

        print('### PARAMETERS ### \n\n')
        offsets = self.parameters['offsets']
        string = '->parameters\n      offsets (len=%d): [' % len(offsets)
        string += ', '.join(offsets)
        string += ']\n'

        if 'date_range' in self.parameters:
            date_range = self.parameters['date_range']
            string += '      date_range: (' + str(date_range[0][0]) + ', ' + str(date_range[0][1]) + ') \n'
        if 'first_hour' in self.parameters:
            first_hour = self.parameters['first_hour']
            string += '      first_hour: (' + str(first_hour[0]) + ', ' + str(first_hour[1]) + ') \n'
        if 'forecast' in self.parameters:
            forecast = self.parameters['forecast']
            if forecast != []:
                string += '      forecast (len=%d): [' % len(forecast)
                for i in forecast:
                    if (len(i) == 2) and isinstance(i, tuple):
                        string += '(%0.2f, %0.2f), ' % i
                    else:
                        string += '--ERROR--, '
                string = string[:-2] + ']\n'
        if 'forecast_d' in self.parameters:
            temp = self.parameters['forecast_d']
            if temp != []:
                string += '      forecast_d (len=%d): [' % len(temp)
                for i in temp:
                    if (len(i) == 2) & (type(i) == tuple):
                        string += '(%0.2f, %0.2f), ' % i
                    else:
                        string = '%s--ERROR--, ' % string
                string = string[:-2] + ']\n'

        if {'type', 'location', 'kind'}.issubset(self.parameters.keys()):
            string += '      type: %s \n      location: %s\n      kind: %s\n\n' % (
            self.parameters['type'], self.parameters['location'], self.parameters['kind'])
        print(string)

        print('\n \n### SECONDARY DATA ###\n\n')
        print("-> dateM (len=%d) ['%s', '%s', '%s', ...]\n" % (
        len(self.dateM), self.dateM[0], self.dateM[1], self.dateM[2]))
        string = 'indexes (len=%d) {' % len(self.indexes)
        incr = 0
        for i in self.indexes:
            if incr >= 6:
                break
            else:
                string += '%d, ' % i
            incr += 1
        string = string[:-2] + '...}\n'
        print(string)

        for i, offset in enumerate(self.parameters['offsets']):
            print('offset: %d' % offset)
            print('-> vectM (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
            len(self.vectM[i]), self.vectM[i][0], self.vectM[i][1], self.vectM[i][2],
            self.vectM[i][3], self.vectM[i][4], self.vectM[i][5]))
            print('-> unif (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
            len(self.unif[i]), self.unif[i][0], self.unif[i][1], self.unif[i][2],
            self.unif[i][3], self.unif[i][4], self.unif[i][5]))
            string = '-> dataM (len=%d) {' % len(self.dataM)
            for key in self.dataM.keys():
                string += key + ', '
            string = string[:-2] + '}\n'
            for key in self.dataM:
                string += '   ' + key + ' (len=%d) [' % len(self.dataM[key][i])
                for j in range(min(6, len(self.dataM[key][i]))):
                    string += '%0.2f, ' % self.dataM[key][i][j]
                string = string[:-2]
                string += '...]\n'
            print(string)
        return None

    # ------------------------------ additional functions ------------------------------------------------

    def plot2D(self):
        if self.dim >= 2:
            plt.figure()
            print(self.vect[:100])
            print(self.vect[:100])
            plt.plot(self.vectM[0], self.vectM[1], '.')

    # this function checks that the parameters have the correct types, length, values...before initialization or update.
    def check_parameters(self, parameters, data):
        self.parameters, self.dim = check_param(parameters, data)
        return None

    # This function creates a 'window' function that tells whether a data point should be considered.
    # This 'window' is then stored in parameters
    def define_window(self):

        if 'window' not in self.parameters:
            dim = self.dim
            parameters = self.parameters
            date_range = None
            if 'date_range' in self.parameters:
                date_range = parameters['date_range']

            def window(l, date, data=None):
                if data is None:
                    data = {}
                if len(l) != dim:
                    raise RuntimeError('l should be of length %d' % dim)

                if any(i is None for i in l):
                    return False

                if 'date_range' in self.parameters:
                    if not any((start_time <= date <= end_time) for start_time, end_time in date_range):
                        return False

                if 'first_hour' in self.parameters:
                    hour = date.hour
                    start_hour, end_hour = parameters['first_hour']
                    if not ((hour - start_hour) % 24 <= (end_hour - start_hour)):
                        return False
                    # b&=(date-3600*(temp[0]-16))%86400<=(temp[1]-temp[0])*3600

                if 'forecast' in self.parameters:
                    for i in range(dim):
                        lower_bound, upper_bound = parameters['forecast'][i]
                        if not (lower_bound <= data['forecast'][i] <= upper_bound):
                            return False

                if 'forecast_d' in self.parameters:
                    for i in range(dim):
                        lower_bound, upper_bound = parameters['forecast_d'][i]
                        if not(lower_bound <= data['forecast_d'][i] <= upper_bound):
                            return False
                return True

            self.parameters['window'] = window
# this function checks that the parameters have the correct types, length, values... before initialization or update.


def check_param(parameters, data):
    # checking mandatory element (offset)

    param_temp = {}

    # print('############ checking parameters ##############')
    if 'offsets' in parameters:
        offsets_fine = all([isinstance(offset, int) for offset in parameters['offsets']])
        dim = len(parameters['offsets'])
        param_temp['offsets'] = sorted(parameters['offsets'])
        if not (offsets_fine and dim > 0):
            raise (RuntimeError('"parameters[\'offsets\']" should be of type \'list of int\' '
                                'and length superior than 0'))
    else:
        raise RuntimeError("'offsets' key in parameters not found ")

    # checking window arguments
    if 'window' in parameters:
        window = parameters['window']
        try:
            if isinstance(window([0.5 for _ in dim], dt.parse('01/01/2000 01:10')), bool):
                param_temp['window'] = parameters['window']
        except TypeError:
            pass

    if 'date_range' in parameters:
        date_range = parameters['date_range']
        try:
            if (isinstance(date_range, (list, tuple)) and len(date_range) == 2
                 and isinstance(date_range[0], str) and isinstance(date_range[1], str)):
                param_temp['date_range'] = [[dt.parse(date) for date in parameters['date_range']]]

            if isinstance(date_range, (list, tuple, set)):
                for date in date_range:
                    if not isinstance(date, (tuple, list)):
                        break
                    if not (len(date) == 2 and isinstance(date[0], str) and isinstance(date[1], str)):
                        break
                else:
                    param_temp['date_range'] = [[dt.parse(start), dt.parse(end)] for start, end in parameters['date_range']]
            else:
                print('Warning: "parameters[\'date_range\']" should be of type (str,str)')
        except ValueError:
            print('Warning: "parameters[\'date_range\']" should be of type (str,str)')

    if 'first_hour' in parameters:
        first_hour = parameters['first_hour']
        if (isinstance(first_hour, tuple)
            and len(first_hour) == 2
            and isinstance(first_hour[0], int) and isinstance(first_hour[1], int)
            and -11 <= first_hour[0] <= first_hour[1] <= 12):

            param_temp['first_hour'] = parameters['first_hour']
        else:
            print('Warning: "parameters[\'first_hour\']" should be of type (int1,int2) with -11<=int1<=int2<=12')

    if 'forecast' in parameters:
        print('\n for: %r \n' % parameters['forecast'])

        if isinstance((parameters['forecast']), list):
            for forecast in parameters['forecast']:
                if not(isinstance(forecast, tuple)
                       and len(forecast) == 2
                       and isinstance(forecast[0], (float, int))
                       and isinstance(forecast[1], (float, int))):
                    break
            else:
                param_temp['forecast'] = parameters['forecast']

        if 'forecast' not in param_temp:
            print('Warning: "parameters[\'forecast\']" should be of type (float,float)')
            print('parameters[\'forecast\']: %r' % parameters['forecast'])
        if 'forecast' in param_temp and ('forecast' not in data):
            del param_temp['forecast']
            print('Warning: you need to give the forecast in the data to use it as a window parameter')

    if 'forecast_d' in parameters:
        if isinstance(parameters['forecast_d'], list):
            for forecast_d in parameters['forecast_d']:
                if not(isinstance(forecast_d, tuple)
                       and len(forecast_d) == 2
                       and isinstance(forecast_d[0], (float, int))
                       and isinstance(forecast_d[1], (float, int))):
                    break
            else:
                param_temp['forecast_d'] = parameters['forecast_d']

        if 'forecast_d' not in param_temp:
            print('Warning: "parameters[\'forecast_d\']" should be of type (float,float)')
        if ('forecast_d' in param_temp) and ('forecast_d' not in data):
            del param_temp['forecast_d']
            print('Warning: you need to give the forecast derivative in the data to use it as a window parameter')

    if 'predicted_day' in parameters:
        if isinstance(parameters['predicted_day'], datetime):
            param_temp['predicted_day'] = parameters['predicted_day']

    # checking informative arguments

    if 'type' in parameters:
        typ = parameters['type']
        if isinstance(typ, str):
            param_temp['type'] = typ
        else:
            print('Warning: "parameters[\'type\']" should be of type string')
    if 'location' in parameters:
        location = parameters['location']
        if isinstance(location, str):
            param_temp['location'] = location
        else:
            print('Warning: "parameters[\'location\']" should be of type string')
    if 'kind' in parameters:
        kind = parameters['kind']
        if isinstance(kind, str):
            param_temp['kind'] = kind
        else:
            print('Warning: "parameters[\'kind\']" should be of type string')

    return param_temp, dim
