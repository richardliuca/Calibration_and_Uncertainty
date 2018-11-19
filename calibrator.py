import os, glob, pyexcel, random, sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
from scipy.signal import butter, medfilt, filtfilt, freqz


def low_pass_filter(x, cutoff = 12.4*10, sampling = 1000, order = 6):
      nyquist = sampling/2
      normalize_freq = (cutoff/nyquist)
      b, a = butter(order, normalize_freq,
                  btype = 'lowpass', analog = False)
      # w, h = freqz(b, a)
      # plt.semilogx(w/(2*np.pi)*sampling, 20*np.log10(abs(h)))
      # plt.axvline(cutoff, color='k')
      # plt.grid(which='both', axis='both')
      # plt.title('Butterworth filter frequency response')
      # plt.xlabel('Frequency [Hz]')
      # plt.ylabel('Amplitude [dB]')
      # plt.show()
      return filtfilt(b, a, x)

def decorator(obj, func , g=lambda x: x, h=lambda y: y):
    def wrapper(*args, **kwargs):
        for file in obj.dataFile:
            xData=np.array(obj.dataBook[file][obj.x], dtype='float64')
            x=g(xData)
            obj.dataBook[file][obj.x]=x.tolist()
            yData=np.array(obj.dataBook[file][obj.y], dtype='float64')
            y=h(yData)
            obj.dataBook[file][obj.y]=y.tolist()
        result=func(*args, **kwargs)
        return result
    return wrapper

# FIXME: Validate initialization in __init__
class Calibrator():
    """
    workDir=working directory
    flag=operation mode (0: all files in workDir, 1: target file)
    target_file=file to calibrate
    x=names for the x-axis column
    y=names for the y-axis column
    z=z or t distrubtion value
    p1=first order polynomial
    p2=second order polynomial
    p3=third order polynomial
    s=sin function
    p=power function
    e=exponential
    l=natural log
    """
    # Reads the data
    # Fitt the data and calculate the model
    # Calculate the Uncertanties in the model fitting

    dataFile=[]
    dataBook={}
    def __init__(self, dir=None, auto=True, target_file=None,
                x=None, y=None, xUnit=None, yUnit=None,
                xUncert=None, yUncert=None, z=1.96, eqn='p1'):
        self.workDir=dir
        self.target=target_file
        _, self.x=self.input_validation(x, type='alpha')
        self.xUnit=xUnit
        _, self.xUncert=self.input_validation(xUncert, type='float')
        _, self.y=self.input_validation(y, type='alpha')
        self.yUnit=yUnit
        _, self.yUncert=self.input_validation(yUncert, type='float')
        self.auto=auto
        _, self.z=self.input_validation(z, type='float')
        _, self.eqn=self.input_validation(eqn, type='alpha')
        if self.auto:
            self.sheetNum=0
        self.request_initialization()
        self.workDir=Path(self.workDir)

    def read_data(self):
        try:
            os.chdir(self.workDir)
        except:
            print('No such directory')
            self.request_initialization()
        fileName=[file for file in os.listdir(self.workDir) if os.path.isfile(file)]
        self.dataFile=[file for file in fileName if 'xlsx' in file or 'csv' in file]
        if not self.target==None and not self.target in self.dataFile:
            return self.error_Handler(0)
        elif not self.target==None and self.target in self.dataFile: # Read target file
            self.dataFile=[self.target]
        elif self.target==None:
            pass
        else:
            return self.error_Handler(0)
        for file in self.dataFile:
            print('Reading {} ...'.format(file))
            if 'xlsx' in file and not self.auto:
                self.sheetNum=self.input_validation(input('Which sheet [0, N]: '), 'int')
            sheet=pyexcel.load(file, sheetname=self.sheetNum, name_columns_by_row=0)
            self.dataBook[file]=sheet.to_dict()
            while not (self.x in self.dataBook[file].keys()) or not (self.y in self.dataBook[file].keys()):
                if not (self.x in self.dataBook[file].keys()):
                    print('X Column not found')
                    self.x=None
                    self.request_initialization()
                if not (self.y in self.dataBook[file].keys()):
                    print('Y column not found')
                    self.y=None
                    self.request_initialization()

    def discrete_fit_data(self, f_custom=None):
        print('Fitting ...')
        f = self.equation_validation(f_custom)
        if self.dataFile==[] or self.dataBook=={}:
            return self.error_Handler(1)
        x=np.array([])
        y=np.array([])
        self.dataBook['xMean']= np.array([])
        self.dataBook['xStd']= np.array([])
        self.dataBook['yUnique']=np.array([])
        self.dataBook['size']=np.array([])
        for file in self.dataFile:
            xFile=np.array(self.dataBook[file][self.x], dtype='float64')
            yFile=np.array(self.dataBook[file][self.y], dtype='float64')
            xMean, xStd, yUnique, size=self.parse_file(x=xFile, y=yFile)
            self.dataBook['xMean']=np.concatenate((self.dataBook['xMean'], xMean))
            self.dataBook['xStd']=np.concatenate((self.dataBook['xStd'], xStd))
            self.dataBook['yUnique']=np.concatenate((self.dataBook['yUnique'],yUnique))
            self.dataBook['size']=np.concatenate((self.dataBook['size'], size))
            x=np.concatenate((x, xFile))
            y=np.concatenate((y, yFile))
        param, covariance=curve_fit(f, self.dataBook['xMean'], self.dataBook['yUnique'], sigma=self.yUncert*np.ones(len(self.dataBook['xMean'])), absolute_sigma=True)
        residuals=y - f(x, *param)
        ss_res=np.sum(residuals**2)
        ss_tot=np.sum((y - np.mean(y))**2)
        r_squared=1 - (ss_res / ss_tot)
        paramStd=np.sqrt(np.diagonal(covariance))
        mean, std, scatterUncert, meanUncert=self.uncertainty_SISO(f, param, self.dataBook['xMean'], self.dataBook['xStd'], self.dataBook['size'])
        message=''
        letter=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        for i in range(len(param)):
            message += letter[i] + '={0:.2f}'.format(param[i]) + '+/- {0:0.2f}'.format(paramStd[i]) + ' '
        print(message)
        print('R^2={0:0.4f}'.format(r_squared))
        self.graph_it(x=x, y=y, f=f, k=param, kStd=self.z*paramStd, msg=message,
                    xMean=self.dataBook['xMean'], yMean=mean,
                    meanUncert=meanUncert, scatterUncert=scatterUncert,
                    xUncert=self.z*self.dataBook['xStd'], yUncert=self.yUncert)
        result = {'param': param,
                'paramStd': paramStd,
                'R2':r_squared,
                'f': f,
                'scatterUncert': scatterUncert.tolist(),
                'meanUncert': meanUncert.tolist()
                }
        return result

    def uncertainty_SISO(self, f=None, k=None, xMean=None, xStd=None, xSize=None):
        print('Simulating ...')
        size=len(xMean)
        outputMean=np.zeros(size, dtype=float)
        outputStd=np.zeros(size, dtype=float)
        for i in range(size):
            progress=(i+1)/size*100
            print('Progress: {0:.02f}%'.format(progress))
            mean=xMean[i]
            std=xStd[i]
            outputMean[i], outputStd[i]=self.monte_carlo_SISO(f, k, mean, std)
        return outputMean, outputStd, self.z*outputStd, self.z*outputStd/xSize

    def monte_carlo_SISO(self, f=None, k=None, mean=None, std=None):
        size=1e5
        random.seed()
        inputSim=[random.gauss(mu=mean, sigma=std) for _ in range(int(size))]
        outputSim=f(np.array(inputSim, dtype='float64'), *k)
        outputMean=np.mean(outputSim, dtype='float64')
        outputStd=np.std(outputSim, dtype='float64')
        return outputMean, outputStd

    def analyzer(self, fx=None, fy=None, inputs={}):
        print('Analyzing ...')
        f = self.equation_validation()
        if self.dataFile==[] or self.dataBook=={}:
            return self.error_Handler(1)
        xStd=inputs['xStd']/self.z
        yStd=inputs['yStd']/self.z
        fxParam=inputs['fxParam']
        fyParam=inputs['fyParam']
        fxUncert=inputs['fxUncert']
        fyUncert=inputs['fyUncert']
        threshold=inputs['threshold']
        result={}
        for file in self.dataFile:
            result.update({file: {}})
            x=self.dataBook[file][self.x]
            y=self.dataBook[file][self.y]
            grouping=self.reduce_transient(x, y, threshold)
            x=sum(grouping['x'], []) # Degree
            y=sum(grouping['y'], []) # Newton
            xxMean, xxStd, xxScatter, xxMeanUncert=self.uncertainty_MISO(fz=fx,
                                                k=fxParam, z=grouping['xMean'],
                                                kStd=fxUncert,
                                                zStd=xStd,
                                                zSize=grouping['xSize'])
            yyMean, yyStd, yyScatter, yyMeanUncert=self.uncertainty_MISO(fz=fy,
                                                k=fyParam,
                                                z=grouping['yMean'],
                                                kStd=fyUncert, zStd=yStd,
                                                zSize=grouping['ySize'])
            progress=(self.dataFile.index(file)+1)/len(self.dataFile)*100
            print('Progress: {0:.02f}%'.format(progress))

            self.graph_it(x=[fx(i, *fxParam) for i in x], y=[fy(i, *fyParam) for i in y],
                        xMean=xxMean, yMean=yyMean,
                        scatterUncert=yyScatter,
                        xUncert=grouping['xStd'] , yUncert=grouping['yStd'])

            result[file].update({'xxMean': xxMean, 'xxStd': xxStd, 'xxScatter': xxScatter, 'xxMeanUncert': xxMeanUncert})
            result[file].update({'yyMean': yyMean, 'yyStd': yyStd, 'yyScatter': yyScatter, 'yyMeanUncert': yyMeanUncert})
        return result

    def uncertainty_MISO(self, fz=None, k=None, z=None, kStd=None, zStd=None, zSize=None):
        outputMean=np.zeros(len(zSize))
        outputStd=np.zeros(len(zSize))
        for i in range(len(zSize)):
            outputMean[i], outputStd[i]=self.monte_carlo_MISO(f=fz, means=[z[i], *k], stds=[zStd, *kStd])
        return outputMean, outputStd, self.z*outputStd, self.z*outputStd/zSize

    def monte_carlo_MISO(self, f=None, means=None, stds=None):
        random.seed()
        size=1E5
        sims=[]
        [sims.append([random.gauss(mu=means[i], sigma=stds[i]) for _ in range(int(size))]) for i in range(len(means))]
        outputSim=list(map(f, *sims))
        outputMean=np.mean(outputSim)
        outputStd=np.std(outputSim)
        return outputMean, outputStd

    def graph_it(self, x=[None], y=[None], f=None, k=[None], kStd=[None], msg=None,
                xMean=[None], yMean=[None], meanUncert=[None], scatterUncert=[None],
                xUncert=None, yUncert=None):
        print('Plotting ...')
        if not(None in x) and not(None in y):
            plt.figure(figsize=(14, 7))
            # plt.plot(x, y, '.',markersize=1, alpha=0.7, label='Experimental Data')
            if not(None in k) and not(None in kStd):
                xInterp=np.arange(min(x), max(x), (max(x) - min(x))/1000)
                plt.plot(xInterp, f(xInterp, *k), 'b-',
                label='Fitted Model' )
                plt.plot(xInterp, f(xInterp, *(k + kStd)), 'r-',
                label='Confidence Band Upper Bound')
                plt.plot(xInterp, f(xInterp, *(k - kStd)), 'r-',
                label='Confidence Band Lower Bound')
                plt.plot(x[0], y[0], 'None', label=msg)
        if not(None in xMean) and not(None in yMean):
            plt.plot(xMean, yMean, 'ro', label='Simulation Mean')
            if not(None in scatterUncert):
                _, capErr, _=plt.errorbar(xMean, yMean, yerr=scatterUncert,
                                        fmt='ko', uplims=True, lolims=True,
                                        label='Scatter Uncertanties')
                capErr[0].set_marker('_')
                capErr[0].set_markersize(20)
                capErr[1].set_marker('_')
                capErr[1].set_markersize(20)
            if  not(None==yUncert):
                _, capUncert, _=plt.errorbar(xMean, yMean, xerr=xUncert, yerr=yUncert,
                                        fmt='o', uplims=True, lolims=True,
                                        label='Measurement Uncertanties')
                capUncert[0].set_marker('_')
                capUncert[0].set_markersize(20)
                capUncert[1].set_marker('_')
                capUncert[1].set_markersize(20)
        plt.xlabel(self.x + ' ' + self.xUnit, size=26)
        plt.ylabel(self.y +  ' ' + self.yUnit, size=26)
        plt.xticks(np.arange(1.1*min(x), 1.1*max(x), (max(x)-min(x))/10))
        plt.yticks(np.arange(1.1*min(y), 1.1*max(y), (max(y)-min(y))/10))
        plt.xticks(size=18)
        plt.yticks(size=18)
        plt.grid(b=True, which='both', axis='both')
        plt.legend(fontsize=16)
        plt.draw()

    def request_initialization(self):
        ''' Initialize required attribute when if not already declared'''
        def is_default(state):
            switcher={'None':False}
            return switcher.get(str(state), True)
        flagDir=is_default(self.workDir)
        flagX=is_default(self.x)
        flagXUncert=is_default(self.xUncert)
        flagY=is_default(self.y)
        flagYUncert=is_default(self.yUncert)
        while not flagX or not flagY or not flagDir  or not flagXUncert or not flagYUncert:
            if not flagDir:
                self.workDir=input('Files directory is: ')
                flagDir=True
            elif not flagX:
                flagX, self.x=self.input_validation(
                                input('Which column is x vector: '), 'str')
                if not flagX:
                    print('X input should be a string')
            elif not flagXUncert:
                flagXUncert, self.xUncert=self.input_validation(
                                input('What is the uncertanties in x: '), 'float')
            elif not flagY:
                flagY, self.y=self.input_validation(
                                input('Which column is y vector: '), 'str')
                if not flagY:
                    print('Y input should be a string')
            elif not flagYUncert:
                flagYUncert, self.yUncert=self.input_validation(
                                input('What is the uncertanties in y: '), 'float')
            else:
                pass

    def input_validation(self, input=None, type=None):
        '''Validates the user inputs with an expected data type'''
        try:
            if input=='' or input==' ':
                return False, None
            elif 'int' in type:
                input=str(input)
                if input.isdigit():
                    return input.isdigit(), int(input)
                else:
                    return False, None
            elif 'float' in type:
                input=str(input)
                digits=input.replace('.', '')
                state, _=self.input_validation(digits, type='int')
                if '.' in input and state:
                    if state:
                        return state, float(input)
                    else:
                        return False, None
                else:
                    return False, None
            elif 'str' in type:
                noSpaces=input.replace(' ', '_')
                noSpecial=noSpaces.replace('_', '')
                if noSpecial.isalpha():
                    return noSpecial.isalpha(), noSpaces
                else:
                    return False, None
            elif 'alpha' in type:
                noUnder=input.replace('_', '')
                noSpaces=noUnder.replace(' ', '')
                noSlash=noSpaces.replace('\\', '')
                noPeriod=noSlash.replace('.', '')
                noColon=noPeriod.replace(':', '')
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

    def equation_validation(self, f_custom=None):
        def model_picker(desire):
            p1=lambda x, a, b: a*x + b
            p2=lambda x, a, b, c: a*x**2.0 + b*x + c
            p3=lambda x, a, b, c, d: a*x**3.0 + b*x**2.0 + c*x +d
            s=lambda x, a, b, c: a*np.sin(b*x + c)
            p=lambda x, a, b: a*x**b
            e=lambda x, a, b, c: a*np.exp(b*x + c)
            l=lambda x, a, b, c: a*np.log(b*x) + c
            models= {'p1': p1, 'p2': p2, 'p3': p3,
                    's': s, 'p': p, 'e': e, 'l': l,
                    'Custom': f_custom, 'None': lambda:'No governing equation found'}
            if f_custom==None:
                return models.get(desire.lower(), lambda:'No governing equation found')
            else:
                return models.get('Custom')
        f=model_picker(str(self.eqn))
        try:
            print(f())
            return self.error_Handler(2)
        except:
            pass
        return f

    def parse_file(self, x=None, y=None):
        yUnique, count=np.unique(y, return_counts=True)
        xMean= np.array([np.mean(x[y==item]) for item in yUnique], dtype='float64')
        xStd=np.array([np.std(x[y==item]) for item in yUnique], dtype='float64')
        return xMean, xStd, yUnique, count

    def reduce_transient(self, x=None, y=None, threshold=None):
        x=np.array(x, dtype='float64')
        y=np.array(y, dtype='float64')
        xReduced=[]
        xMean=[]
        xStd=[]
        xSize=[]
        yReduced=[]
        yMean=[]
        yStd=[]
        ySize=[]
        binsFiltered=[]
        occur, bins=np.histogram(x, bins='fd', density=False)
        [binsFiltered.append([bins[i], bins[i+1]]) for i in range(len(bins)-1) if occur[i] >=threshold]
        for i in range(len(binsFiltered)):
            lowerBound=binsFiltered[i][0]
            upperBound=binsFiltered[i][1]
            xEntries=x[np.logical_and(x>=lowerBound, x<=upperBound)].tolist()
            yEntries=y[np.logical_and(x>=lowerBound, x<=upperBound)].tolist()
            xReduced.append(xEntries)
            xMean.append(np.mean(xEntries))
            xStd.append(np.std(xEntries))
            xSize.append(len(xEntries))
            yReduced.append(yEntries)
            yMean.append(np.mean(yEntries))
            yStd.append(np.std(yEntries))
            ySize.append(len(yEntries))
        return {'x':xReduced, 'xMean':xMean, 'xStd':xStd, 'xSize':xSize,
                'y':yReduced, 'yMean':yMean, 'yStd':yStd, 'ySize':ySize}

    def error_Handler(self, code=None):
        if code==0:
            print('File not found')
        elif code==1:
            print('Read data first by calling read_data')
        elif code==2:
            print('No data was fitted')
        elif code==3:
            print('Validation failed')
        elif code==4:
            print('Nothing to plot')
