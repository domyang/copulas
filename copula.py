import matplotlib.pyplot as plt
# from pyemd import emd
from pylab import rcParams
import numpy as np
import scipy as sc
import dateutil.parser as dt
from dateutil.relativedelta import relativedelta
import spline
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
    vect = []  # a list of float: the real data to be considered
    date = []  # a list (of ints) of ordered date
    data = {}  # a dictionary of lists, additional parameters considered during the creation of the window

    # dictionary defining 3 things:
    #                           the offsets (list)
    #                           the window (numerous arguments)
    #                           the kind of data that has been given in argument
    parameters = {"offsets": [], 'predicted_day': '2000-01-01 00:00',
                  'date_range': ('', ''), 'first_hour': (-11, 12), 'forecast': [], 'forecast_d': [],
                  'type': 'Wind', 'location': 'total', 'kind': 'error'}

    # dictionary saying which parameters to take into account
    param_actual = {'window': False, 'date_range': False, 'first_hour': False, 'forecast': False, 'forecast_d': False}

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
    def __init__(self, arg_date, vect, data, parameters, var_function=None):

        if type(arg_date[0]) is int:
            date = arg_date
        else:
            date = [int(dt.parse(i).timestamp()) for i in arg_date]
        data_keys = data.keys()

        print('preparing the data')
        if len(date) == len(vect):

            vect_temp = []
            date_temp = []
            data_temp = {}
            for i in data_keys:
                data_temp[i] = []
            dt_object = dt.parse(arg_date[0])

            vect_it = iter(vect)
            date_it = iter(date)
            data_it = {}
            for key in data.keys():
                data_it[key] = iter(data[key])

            cur = date[0]
            while (cur < date[-1]):
                index = next(date_it)

                while cur < index:
                    vect_temp.append(None)
                    for key in data.keys():
                        data_temp[key].append(None)
                    date_temp.append(cur)
                    cur = cur + 3600

                while cur > index:
                    print(
                        'Warning: possibly twice the same entry at date %s' % dt_object.fromtimestamp(index).__str__())
                    index = next(date_it)
                    next(vect_it)
                    for key in data.keys():
                        next(data_it[key])

                if cur == index:
                    date_temp.append(index)
                    vect_temp.append(next(vect_it))
                    for key in data.keys():
                        data_temp[key].append(next(data_it[key]))

                cur = cur + 3600

            self.vect = vect_temp
            self.date = date_temp
            self.length = len(vect_temp)
            self.data = data_temp
            self.var_function = var_function

        else:
            raise (
            RuntimeError('vect and date should be same length, (here %d,%d respectively)' % (self.length, len(date))))

        self.update(parameters)

        print('copula initialized')

    ### the update functions update the window to match the parameters arguments
    # returns nothing
    # arguments:
    #   () parameters: a dictionary specifying the window
    #   () data: additional data to be taken into account
    def update(self, parameters):

        print('checking the parameters arguments')
        data_temp = self.data
        self.check_parameters(parameters, data_temp)
        self.define_window()

        # initializing the temporary variables
        dim = self.dim
        vectM_temp = [[] for i in range(dim)]
        indexes = []
        dataM_temp = {}
        dateM_temp = []

        data_keys = self.data.keys()

        for key in data_keys:
            dataM_temp[key] = [[] for i in range(dim)]
        max_offset = max(parameters['offsets']) + 1

        # buffers for the next loop
        buff = [0 for i in range(max_offset)]
        buff_ind = 0
        data_buff = {}

        for i in data_keys:
            data_buff[i] = buff.copy()

        # variables needed
        vect_iter = iter(self.vect)
        date_iter = iter(self.date)
        data_iter = {}
        for key in data_keys:
            data_iter[key] = iter(self.data[key])

        print('selecting the desired points')
        # selecting the desired points

        for ind in range(self.length):
            buff[buff_ind] = next(vect_iter)
            for key in data_keys:
                data_buff[key][buff_ind] = next(data_iter[key])
            buff_ind = (buff_ind + 1) % max_offset

            if ind < max_offset - 1:
                continue

            cur = next(date_iter)

            temp = [[], {}]
            for i in parameters['offsets']:
                temp[0].append(buff[(i + buff_ind) % max_offset])

            for key in data_keys:
                temp[1][key] = []
                for i in parameters['offsets']:
                    temp[1][key].append(data_buff[key][(i + buff_ind) % max_offset])

            if self.parameters['window'](temp[0], cur, temp[1]):

                for i in range(dim):
                    vectM_temp[i].append(temp[0][i])
                indexes.append(ind - max_offset)
                for key in temp[1].keys():
                    for i in range(dim):
                        dataM_temp[key][i].append(temp[1][key][i])
                dateM_temp.append(cur)

        self.dateM = dateM_temp
        self.dataM = dataM_temp
        self.lengthM = len(vectM_temp[0])

        # if this is errors in Solar power forecasts, applying reverse variance transformation to get the correct values in vectM
        # (variances were scaled to 1 for each solar hour, now we multiply by the variance corresponding
        # to the solar hour of 'predicted_day' with offset)
        if (parameters['type'] == 'Solar') & (parameters['kind'] == 'error'):
            if ('hour_sol' in data_keys) & (self.var_function is not None):
                if not 'predicted_day' in parameters.keys():

                    if 'date_range' in parameters.keys():
                        parameters['predicted_day'] = parameters['date_range'][1]
                    else:
                        parameters['predicted_day'] = '2000-01-01 00:00'

                time = dt.parse(parameters['predicted_day'])
                var_sol_temp = []
                rise, sete = ut.sun_rise_set('%d-%d-%d' % (time.year, time.month, time.day))
                for i in parameters['offsets']:
                    var_sol_temp.append(self.var_function(
                        ((i + time.hour + time.minute / 60 + time.second / 3600) - rise) / (sete - rise) * 12))
                for i in range(dim):
                    def f(x):
                        return x * var_sol_temp[i]

                    vectM_temp[i] = list(map(f, vectM_temp[i]))

        self.vectM = vectM_temp
        # creating the copula variables
        print('creating the copula points')
        indexes = set(indexes)

        self.indexes = indexes
        self.unif = ut.uniforms(vectM_temp, rand=False)

    ### pprint doesn't plot: it just list all the attributes:
    def pprint(self):

        print('\n### SCALAR PARAMETERS ### \n \n length= %d , dim= %d , lengthM= %d \n\n' % (
        self.length, self.dim, self.lengthM))
        print('### PRIMARY DATA ### \n \n-> vect (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
        len(self.vect), self.vect[0], self.vect[1], self.vect[2], self.vect[3], self.vect[4], self.vect[5]))
        print('-> date (len=%d) [\'%s\', \'%s\', \'%s\', ...]\n' % (
        len(self.date), self.date[0], self.date[1], self.date[2]))
        string = '-> data (len=%d) {' % len(self.data)
        for key in self.data.keys():
            string += '\'%s\', ' % key
            string = string[:-2]
        string += '}\n'
        for key in self.data.keys():
            string += '   ' + key + ' (len=%d) [' % len(self.data[key])
            for i in range(min(6, len(self.data[key]))):
                string += '%0.2f, ' % self.data[key][i]
            string += '...]\n'
        print(string)

        print('### PARAMETERS ### \n\n')
        temp = self.parameters['offsets']
        string = '->parameters\n      offsets (len=%d): [' % len(temp)
        for i in temp:
            string += '%d, ' % i
        string = string[:-2]
        string += ']\n'

        if 'date_range' in self.parameters.keys():
            temp = self.parameters['date_range']
            string += '      date_range: (' + str(temp[0][0]) + ', ' + str(temp[0][1]) + ') \n'
        if 'first_hour' in self.parameters.keys():
            temp = self.parameters['first_hour']
            string += '      first_hour: (' + str(temp[0]) + ', ' + str(temp[1]) + ') \n'
        if 'forecast' in self.parameters.keys():
            temp = self.parameters['forecast']
            if temp != []:
                string += '      forecast (len=%d): [' % len(temp)
                for i in temp:
                    if (len(i) == 2) & (type(i) == tuple):
                        string += '(%0.2f, %0.2f), ' % i
                    else:
                        string += '--ERROR--, '
                string = string[:-2] + ']\n'
        if 'forecast_d' in self.parameters.keys():
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
        print('param_actual')
        print('      %r' % self.param_actual)

        print('\n \n### SECONDARY DATA ###\n\n')
        print('-> dateM (len=%d) [\'%s\', \'%s\', \'%s\', ...]\n' % (
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

        for offset in range(len(self.parameters['offsets'])):
            print('offset: %d' % self.parameters['offsets'][offset])
            print('-> vectM (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
            len(self.vectM[offset]), self.vectM[offset][0], self.vectM[offset][1], self.vectM[offset][2],
            self.vectM[offset][3], self.vectM[offset][4], self.vectM[offset][5]))
            print('-> unif (len=%d) [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f...]\n' % (
            len(self.unif[offset]), self.unif[offset][0], self.unif[offset][1], self.unif[offset][2],
            self.unif[offset][3], self.unif[offset][4], self.unif[offset][5]))
            string = '-> dataM (len=%d) {' % len(self.dataM)
            for key in self.dataM.keys():
                string += key + ', '
            string = string[:-2] + '}\n'
            for key in self.dataM.keys():
                string += '   ' + key + ' (len=%d) [' % len(self.dataM[key][offset])
                for i in range(min(6, len(self.dataM[key][offset]))):
                    string += '%0.2f, ' % self.dataM[key][offset][i]
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

    ### this function checks that the parameters have the correct types, length, values... before initialization or update.
    def check_parameters(self, parameters, data):
        self.parameters, self.param_actual, self.dim = check_param(parameters, data)
        return None

    ### This function creates a 'window' function that tells whether a data point should be considered.
    ### This 'window' is then stored in parameters
    def define_window(self):

        dt_obj = dt.parse('2000-1-1 00:00')

        if not self.param_actual['window']:
            dim = self.dim
            parameters = self.parameters
            param_actual = self.param_actual
            date_range = None
            if param_actual['date_range']:
                date_range = [[dt.parse(i).timestamp() for i in j] for j in parameters['date_range']]

            def f(l, date, data={}):
                b = True
                if len(l) == dim:
                    no_none = True;
                    for i in l:
                        if i is None:
                            no_none = False
                            break
                    if no_none:

                        if param_actual['date_range']:
                            good_date = False
                            for date_tp in date_range:
                                good_date |= (date_tp[0] <= date <= date_tp[1])
                                if good_date:
                                    break
                            b &= good_date

                        if param_actual['first_hour']:
                            hour = dt_obj.fromtimestamp(date).hour
                            temp = parameters['first_hour']
                            b &= (hour - temp[0]) % 24 <= (temp[1] - temp[0])
                            # b&=(date-3600*(temp[0]-16))%86400<=(temp[1]-temp[0])*3600

                        if param_actual['forecast']:
                            for i in range(dim):
                                b &= parameters['forecast'][i][0] <= data['forecast'][i] <= parameters['forecast'][i][1]

                        if param_actual['forecast_d']:
                            for i in range(dim):
                                b &= parameters['forecast_d'][i][0] <= data['forecast_d'][i] <= \
                                     parameters['forecast_d'][i][1]
                    if not no_none:
                        return False
                else:
                    raise (RuntimeError('l should be of length %d' % dim))
                return b

            self.parameters['window'] = f
            self.param_actual['window'] = True


### this function checks that the parameters have the correct types, length, values... before initialization or update.
def check_param(parameters, data):
    keys = parameters.keys()

    # checking mandatory element (offset)

    param_temp = {}
    param_actual = {'window': False, 'date_range': False, 'first_hour': False, 'kind': False, 'forecast': False,
                    'forecast_d': False}
    dim = 0

    print('############ checking parameters ##############')
    if 'offsets' in keys:
        offsets_fine = ut.is_list_of(parameters['offsets'], child=int)[0]
        if offsets_fine:
            dim = len(parameters['offsets'])
            offsets_fine &= dim > 0
            parameters['offsets'].sort()
            param_temp['offsets'] = parameters['offsets']
        if not offsets_fine:
            raise (RuntimeError('"parameters[\'offsets\']" should be of type \'list of int\' '
                                'and length superior than 0'))
    else:
        raise (RuntimeError('\'offsets\' key in parameters not found '))

    # checking window arguments


    if 'window' in keys:
        temp = parameters['window']
        try:
            if type(temp([0.5 for i in dim], '01/01/2000 01:10')) == bool:
                param_actual['window'] = True
                param_temp['window'] = parameters['window']
        except TypeError:
            pass

    if 'date_range' in keys:
        temp = parameters['date_range']
        try:
            if type(temp) in {list, tuple}:
                if len(temp) == 2:
                    if type(temp[0]) == type(temp[1]) == str:
                        param_actual['date_range'] = True
                        param_temp['date_range'] = [parameters['date_range']]
            if type(temp) in {list, tuple, set}:
                b = True
                for i in temp:
                    if not type(i) in {tuple, list}:
                        b = False
                    correct = False
                    if len(i) == 2:
                        if type(i[0]) == type(i[1]) == str:
                            correct = True
                    b &= correct
                    if b == False:
                        break
                if b:
                    param_actual['date_range'] = True
                    param_temp['date_range'] = parameters['date_range']

            if not param_actual['date_range']:
                print('Warning: "parameters[\'date_range\']" should be of type (str,str)')
        except ValueError:
            print('Warning: "parameters[\'date_range\']" should be of type (str,str)')
    if 'first_hour' in keys:
        temp = parameters['first_hour']
        if type(temp) == tuple:
            if len(temp) == 2:
                if type(temp[0]) == type(temp[1]) == int:
                    if -11 <= temp[0] <= temp[1] <= 12:
                        param_actual['first_hour'] = True
                        param_temp['first_hour'] = parameters['first_hour']
        if not param_actual['first_hour']:
            print('Warning: "parameters[\'first_hour\']" should be of type (int1,int2) with -11<=int1<=int2<=12')
    if 'forecast' in keys:

        print('\n for: %r \n' % parameters['forecast'])

        b = type(parameters['forecast']) == list
        if b:
            for temp in parameters['forecast']:
                b_bis = False
                if type(temp) == tuple:
                    if len(temp) == 2:
                        if (type(temp[0]) in {float, int}) & (type(temp[1]) in {float, int}):
                            b_bis = True
                b &= b_bis
        param_actual['forecast'] = b
        if b:
            param_temp['forecast'] = parameters['forecast']
        else:
            print('Warning: "parameters[\'forecast\']" should be of type (float,float)')
            print('parameters[\'forecast\']: %r' % parameters['forecast'])
        if param_actual['forecast'] & (not 'forecast' in data.keys()):
            param_actual['forecast'] = False
            print('Warning: you need to give the forecast in the data to use it as a window parameter')
    if 'forecast_d' in keys:
        b = type(parameters['forecast_d']) == list
        if b:
            for temp in parameters['forecast_d']:
                b_bis = False
                if type(temp) == tuple:
                    if len(temp) == 2:
                        if (type(temp[0]) in {float, int}) & (type(temp[1]) in {float, int}):
                            b_bis = True
                b &= b_bis
        param_actual['forecast_d'] = b
        if b:
            param_temp['forecast_d'] = parameters['forecast_d']
        else:
            print('Warning: "parameters[\'forecast_d\']" should be of type (float,float)')
        if param_actual['forecast_d'] & (not 'forecast_d' in data.keys()):
            param_actual['forecast_d'] = False
            print('Warning: you need to give the forecast derivative in the data to use it as a window parameter')

    # checking informative arguments

    if 'type' in keys:
        temp = parameters['type']
        if type(temp) == str:
            param_temp['type'] = parameters['type']
        else:
            print('Warning: "parameters[\'type\']" should be of type string')
    if 'location' in keys:
        temp = parameters['location']
        if type(temp) == str:
            param_temp['location'] = parameters['location']
        else:
            print('Warning: "parameters[\'location\']" should be of type string')
    if 'kind' in keys:
        temp = parameters['kind']
        if type(temp) == str:
            param_temp['kind'] = parameters['kind']
        else:
            print('Warning: "parameters[\'kind\']" should be of type string')

    return (param_temp, param_actual, dim)
