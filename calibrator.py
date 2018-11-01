import os, glob, pyexcel, random
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path


# FIXME: Add x y data modifer (decorator)
# FIXME: Validate initialization in __init__
class calibrator():
    """
    workDir = working directory
    flag = operation mode (0: all files in workDir, 1: target file)
    target_file = file to calibrate
    x = names for the x-axis column
    y = names for the y-axis column
    z = z or t distrubtion value
    p1 = first order polynomial
    p2 = second order polynomial
    p3 = third order polynomial
    s = sin function
    p = power function
    e = exponential
    l = natural log
    """
    # Reads the data
    # Fitt the data and calculate the model
    # Calculate the Uncertanties in the model fitting
    dataFile = []
    dataBook = {}
    def __init__(self, dir = None, auto = True, target_file = None,
                x = None, y = None, xUnit = None, yUnit = None,
                xUncert = None, yUncert = None,
                z = 1.96, eqn = 'p1'):
        self.workDir = dir
        self.target = target_file
        _, self.x = self.input_validation(x, type = 'alpha')
        self.xUnit = xUnit
        _, self.xUncert = self.input_validation(xUncert, type = 'float')
        _, self.y = self.input_validation(y, type = 'alpha')
        self.yUnit = yUnit
        _, self.yUncert = self.input_validation(yUncert, type = 'float')
        self.auto = auto
        _, self.z = self.input_validation(z, type = 'float')
        _, self.eqn = self.input_validation(eqn, type = 'alpha')
        if self.auto == True:
            self.sheetNum = 0;
        self.request_initialization()
        self.workDir = Path(self.workDir)

    def read_data(self):
        try:
            os.chdir(self.workDir)
        except:
            print('No such directory')
            self.request_initialization()
        fileName = [file for file in os.listdir(self.workDir) if os.path.isfile(file)]
        self.dataFile = [file for file in fileName if 'xlsx' in file or 'csv' in file]
        if not self.target == None and not self.target in self.dataFile:
            return self.error_Handler(0)
        elif not self.target == None and self.target in self.dataFile: # Read target file
            self.dataFile = [self.target]
        elif self.target == None:
            pass
        else:
            self.error_Handler(0)
        for file in self.dataFile:
            print('Reading {} ...'.format(file))
            if 'xlsx' in file and self.auto == False:
                self.sheetNum = self.input_validation(input('Which sheet [0, N]: '), 'int')
            sheet = pyexcel.load(file, sheetname = self.sheetNum, name_columns_by_row = 0)
            self.dataBook[file] = sheet.to_dict()
            while not (self.x in self.dataBook[file].keys()) or not (self.y in self.dataBook[file].keys()):
                if not (self.x in self.dataBook[file].keys()):
                    print('X Column not found')
                    self.x = None
                    self.request_initialization()
                if not (self.y in self.dataBook[file].keys()):
                    print('Y column not found')
                    self.y = None
                    self.request_initialization()

    def discrete_fit_data(self, f_custom  = None):
        print('Fitting ...')
        def model_picker(desire):
            p1 = lambda x, a, b: a*x + b
            p2 = lambda x, a, b, c: a*x**2.0 + b*x + c
            p3 = lambda x, a, b, c, d: a*x**3.0 + b*x**2.0 + c*x +d
            s = lambda x, a, b, c: a*np.sin(b*x + c)
            p = lambda x, a, b: a*x**b
            e = lambda x, a, b, c: a*np.exp(b*x + c)
            l = lambda x, a, b, c: a*np.log(b*x) + c
            models= {'p1': p1, 'p2': p2, 'p3': p3,
                    's': s, 'p': p, 'e': e, 'l': l,
                    'Custom': f_custom, 'None': lambda:'No governing equation found'}
            if f_custom == None:
                return models.get(desire.lower(), lambda:'No governing equation found')
            else:
                return models.get('Custom')
        f = model_picker(str(self.eqn))
        try:
            print(f())
            return self.error_Handler(2)
        except:
            pass
        if self.dataFile == [] or self.dataBook == {}:
            return self.error_Handler(1)
        x = np.array([])
        y = np.array([])
        self.dataBook['xMean'] =  np.array([])
        self.dataBook['xStd'] =  np.array([])
        self.dataBook['yUnique'] = np.array([])
        self.dataBook['size'] = np.array([])
        for file in self.dataFile:
            xFile = np.array(self.dataBook[file][self.x], dtype = float)
            yFile = np.array(self.dataBook[file][self.y], dtype = float)
            xMean, xStd, yUnique, size = self.parse_file(x = xFile, y = yFile)
            self.dataBook['xMean'] = np.concatenate((self.dataBook['xMean'], xMean))
            self.dataBook['xStd'] = np.concatenate((self.dataBook['xStd'], xStd))
            self.dataBook['yUnique'] = np.concatenate((self.dataBook['yUnique'],yUnique))
            self.dataBook['size'] = np.concatenate((self.dataBook['size'], size))
            x = np.concatenate((x, xFile))
            y = np.concatenate((y, yFile))
        param, covariance = curve_fit(f, self.dataBook['xMean'], self.dataBook['yUnique'], sigma = self.yUncert*np.ones(len(self.dataBook['xMean'])), absolute_sigma = True)
        residuals = y - f(x, *param)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        paramStd = np.sqrt(np.diagonal(covariance))
        mean, std, scatterUncert, meanUncert = self.uncertainty_SISO(f, param, self.dataBook['xMean'], self.dataBook['xStd'], self.dataBook['size'])
        message = ''
        letter = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        for i in range(len(param)):
            message += letter[i] + ' = {0:.2f}'.format(param[i]) + '+/- {0:0.2f}'.format(paramStd[i]) + ' '
        print(message)
        print('R^2 = {0:0.4f}'.format(r_squared))
        self.graph_it(x = x, y = y, f = f, k = param, kStd = self.z*paramStd, msg = message,
                    xMean = self.dataBook['xMean'], yMean = mean,
                    meanUncert = meanUncert, scatterUncert = scatterUncert,
                    xUncert = self.xUncert, yUncert = self.yUncert)

    def uncertainty_SISO(self, f = None, k = None, xMean = None, xStd = None, xSize = None):
        print('Simulating ...')
        size = len(xMean)
        outputMean = np.zeros(size, dtype=float)
        outputStd = np.zeros(size, dtype = float)
        for i in range(size):
            mean = xMean[i]
            std = xStd[i]
            outputMean[i], outputStd[i] = self.monte_carlo_SISO(f, k, mean, std)
        return outputMean, outputStd, self.z*outputStd, self.z*outputStd/xSize

    def monte_carlo_SISO(self, f = None, k = None, mean = None, std = None):
        size = 1e6
        random.seed()
        inputSim = [random.gauss(mu = mean, sigma = std) for _ in range(int(size))]
        outputSim = f(np.array(inputSim), *k)
        outputMean = np.mean(outputSim, dtype = np.float64)
        outputStd = np.std(outputSim, dtype = np.float64)
        return outputMean, outputStd

    def graph_it(self, x = None, y = None, f = None, k = None, kStd = None, msg = None,
                xMean = None, yMean = None, meanUncert = None, scatterUncert = None,
                xUncert = None, yUncert = None):
        print('Plotting ...')
        plt.figure(figsize = (14, 7))
        plt.plot(x, y, '.',markersize = 4, label = 'Experimental Data')

        xInterp = np.arange(min(x), max(x), (max(x) - min(x))/1000)
        plt.plot(xInterp, f(xInterp, *k), 'b-',
        label = 'Fitted Model' )
        plt.plot(xInterp, f(xInterp, *(k + kStd)), 'r-',
        label = 'Confidence Band Upper Bound')
        plt.plot(xInterp, f(xInterp, *(k - kStd)), 'r-',
        label = 'Confidence Band Lower Bound')
        plt.plot(x[0], y[0], 'None', label = msg)

        plt.plot(xMean, yMean, 'ro', label = 'Simulation Mean')
        plt.errorbar(xMean, yMean, yerr = scatterUncert, fmt = 'k*',
        uplims = True, lolims = True, label = 'Scatter Uncertanties')

        plt.errorbar(xMean, yMean, xerr = xUncert, yerr = yUncert,
        fmt = '*', uplims = True, lolims = True, label = 'Measurement Uncertanties')

        plt.xlabel(self.x + ' ' + self.xUnit)
        plt.ylabel(self.y +  ' ' + self.yUnit)
        plt.grid(b = True, which = 'both')
        plt.legend(fontsize = 12)
        plt.draw()

    def request_initialization(self):
        ''' Initialize required attribute when if not already declared'''
        def is_default(state):
            switcher = {'None':False}
            return switcher.get(str(state), True)
        flagDir = is_default(self.workDir)
        flagX = is_default(self.x)
        flagXUncert = is_default(self.xUncert)
        flagY = is_default(self.y)
        flagYUncert = is_default(self.yUncert)
        while not flagX or not flagY or not flagDir  or not flagXUncert or not flagYUncert:
            if not flagDir:
                self.workDir = input('Files directory is: ')
                flagDir = True
            elif not flagX:
                flagX, self.x = self.input_validation(
                                input('Which column is x vector: '), 'str')
                if not flagX:
                    print('X input should be a string')
            elif not flagXUncert:
                flagXUncert, self.xUncert = self.input_validation(
                                input('What is the uncertanties in x: '), 'float')
            elif not flagY:
                flagY, self.y = self.input_validation(
                                input('Which column is y vector: '), 'str')
                if not flagY:
                    print('Y input should be a string')
            elif not flagYUncert:
                flagYUncert, self.yUncert = self.input_validation(
                                input('What is the uncertanties in y: '), 'float')
            else:
                pass

    def input_validation(self, input = None, type = None):
        '''Validates the user inputs with an expected data type'''
        try:
            if input == '' or input == ' ':
                return False, None
            elif 'int' in type:
                input = str(input)
                if input.isdigit():
                    return input.isdigit(), int(input)
                else:
                    return False, None
            elif 'float' in type:
                input = str(input)
                digits = input.replace('.', '')
                state, _ = self.input_validation(digits, type = 'int')
                if '.' in input and state:
                    if state:
                        return state, float(input)
                    else:
                        return False, None
                else:
                    return False, None
            elif 'str' in type:
                noSpaces = input.replace(' ', '_')
                noSpecial = noSpaces.replace('_', '')
                if noSpecial.isalpha():
                    return noSpecial.isalpha(), noSpaces
                else:
                    return False, None
            elif 'alpha' in type:
                noUnder = input.replace('_', '')
                noSpaces = noUnder.replace(' ', '')
                noSlash = noSpaces.replace('\\', '')
                noPeriod = noSlash.replace('.', '')
                noColon = noPeriod.replace(':', '')
                if noColon.isalnum():
                    return noColon.isalnum(), input
                else:
                    return False, None
            else:
                self.error_Handler(3)
                return False, None
        except:
            print('Data conversion Failed, Input not valid')
            return False, None

    def parse_file(self, x = None, y = None):
        yUnique, count = np.unique(y, return_counts = True)
        xMean =  np.array([np.mean(x[y == item]) for item in yUnique])
        xStd = np.array([np.std(x[y == item]) for item in yUnique])
        return xMean, xStd, yUnique, count


    def error_Handler(self, code = None):
        if code == 0:
            print('File not found')
        elif code == 1:
            print('Read data first by calling read_data')
        elif code == 2:
            print('No data was fitted')
        elif code == 3:
            print('Validation failed')
        elif code == 4:
            print('Nothing to plot')

    def debug(self, msg = None):
        print('Working in {}'.format(self.workDir))
        print('Working in {}'.format(self.dataBook[self.target]))


if __name__ == '__main__':
    workingDirectory = r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\10.26.2018'
    settings = {'dir': workingDirectory,
                'x': 'Ch1_Fz', 'y': 'Mass',
                'xUnit': '[volts]', 'yUnit': '[gram]',
                'xUncert': 0.0, 'yUncert': 0.2,
                'target_file': 'strain_gauge_cal.csv', 'auto': True,
                'eqn': 'p1'}
    customEquation = None
    bot1 = calibrator(**settings)
    bot1.read_data()
    bot1.discrete_fit_data(customEquation)
    workingDirectory = r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\10.26.2018\Angle Attack Cali'
    settings = {'dir': workingDirectory,
                'x': 'Ch4_a_control', 'y': 'Angle_of_Attack',
                'xUnit': '[volts]', 'yUnit': '[degree]',
                'xUncert': 0.0, 'yUncert': 0.1,
                'auto': True,
                'eqn': 'p1'}
    customEquation = None
    bot2 = calibrator(**settings)
    bot2.read_data()
    bot2.discrete_fit_data(customEquation)

    plt.show()
    plt.close()
