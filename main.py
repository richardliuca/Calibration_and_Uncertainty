from calibrator import Calibrator, decorator
import matplotlib.pyplot as plt
import numpy as np
import os, collections, pyexcel

armLength=14/100
maxWeight=150/1000*9.81
momentArmTilt=4.2*np.pi/180

armLengthUncert=0.5/100
weightUncert=1/1000*9.81
levelUncert=0.1

dragBias=3.550587587-0.911987642

savingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\Result'
os.chdir(savingDirectory)
files=[file for file in os.listdir(savingDirectory) if os.path.isfile(file)]
if 'results.xlsx' in files:
    SG=pyexcel.load('results.xlsx', sheetname='SG', name_columns_by_row=0)
    resultSG=dict(SG.to_dict())
    resultSG.update({'f':lambda x, a, b: a*x + b ,
                    'param':resultSG['Parameter'],
                    'paramStd':resultSG['Parameter_Std'],
                    'scatterUncert':[resultSG['Max95%Scatter_Uncertainty'][0]],
                    'meanUncert':[resultSG['Max95%Mean_Uncertainty'][0]]})

    AOA=pyexcel.load('results.xlsx', sheetname='AOA', name_columns_by_row=0)
    resultAOA=dict(AOA.to_dict())
    resultAOA.update({'f':lambda x, a, b: a*x + b,
                    'param':resultAOA['Parameter'],
                    'paramStd':resultAOA['Parameter_Std'],
                    'scatterUncert':[resultAOA['Max95%Scatter_Uncertainty'][0]],
                    'meanUncert':[resultAOA['Max95%Mean_Uncertainty'][0]]})

    Drag=pyexcel.load('results.xlsx', sheetname='Drag', name_columns_by_row=0)
    resultDrag=dict(Drag.to_dict())
    resultDrag.update({'f':lambda x, a, b: a*x + b,
                    'param':resultDrag['Parameter'],
                    'paramStd':resultDrag['Parameter_Std'],
                    'scatterUncert':[resultDrag['Max95%Scatter_Uncertainty'][0]],
                    'meanUncert':[resultDrag['Max95%Mean_Uncertainty'][0]]})

    Moment=pyexcel.load('results.xlsx', sheetname='MomentTotal', name_columns_by_row=0)
    resultMoment=dict(Moment.to_dict())
    resultMoment.update({'f':lambda x, a, b: a*x + b,
                    'param':resultMoment['Parameter'],
                    'paramStd':resultMoment['Parameter_Std'],
                    'scatterUncert':[resultMoment['Max95%Scatter_Uncertainty'][0]],
                    'meanUncert':[resultMoment['Max95%Mean_Uncertainty'][0]]})

else:
    # NOTE: SG Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\10.26.2018'
    settings={'dir': workingDirectory,
                'x': 'Ch1_Fz', 'y': 'Weight',
                'xUnit': '[volts]', 'yUnit': '[Newton]',
                'xUncert': 0.0, 'yUncert': weightUncert,
                'target_file': 'strain_gauge_cal.csv', 'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botSG=Calibrator(**settings)
    botSG.read_data()
    discrete_fit_data_SG=decorator(botSG, botSG.discrete_fit_data, h=lambda y: -y/1000*9.81)
    resultSG=discrete_fit_data_SG(customEquation)
    SGData=[
            ['Parameter', 'Parameter_Std', 'R2', 'Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultSG['param'][0], resultSG['paramStd'][0], resultSG['R2'], max(resultSG['scatterUncert']), max(resultSG['meanUncert'])],
            [resultSG['param'][1], resultSG['paramStd'][1], '', '', '']
            ]

    # NOTE: AOA Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\10.26.2018\Angle Attack Cali'
    settings={'dir': workingDirectory,
                'x': 'Ch4_a_control', 'y': 'Angle_of_Attack',
                'xUnit': '[volts]', 'yUnit': '[degree]',
                'xUncert': 0.0, 'yUncert': levelUncert,
                'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botAOA=Calibrator(**settings)
    botAOA.read_data()
    discrete_fit_data_AOA=decorator(botAOA, botAOA.discrete_fit_data)
    resultAOA=discrete_fit_data_AOA(customEquation)
    AOAData=[
            ['Parameter', 'Parameter_Std', 'R2','Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultAOA['param'][0], resultAOA['paramStd'][0], resultAOA['R2'], max(resultAOA['scatterUncert']), max(resultAOA['meanUncert'])],
            [resultAOA['param'][1], resultAOA['paramStd'][1], '', '', '']
            ]

    # NOTE: Drag Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\WS16_VWT-V10-01_calibration\Fx-drag'
    settings={'dir': workingDirectory,
                'x': 'Ch0_Fx', 'y': 'Drag',
                'xUnit': '[volts]', 'yUnit': '[Newton]',
                'xUncert': 0.0, 'yUncert': weightUncert,
                'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botDrag=Calibrator(**settings)
    botDrag.read_data()
    bias=3.550587587-0.911987642
    discrete_fit_data_Drag=decorator(botDrag, botDrag.discrete_fit_data, g=lambda x: x-dragBias, h=lambda y: y/1000*9.81)
    resultDrag=discrete_fit_data_Drag(customEquation)
    DragData=[
            ['Parameter', 'Parameter_Std', 'R2', 'Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultDrag['param'][0], resultDrag['paramStd'][0], resultDrag['R2'], max(resultDrag['scatterUncert']), max(resultDrag['meanUncert'])],
            [resultDrag['param'][1], resultDrag['paramStd'][1], '', '', '']
            ]

    # NOTE: CCW Moment Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\WS16_VWT-V10-01_calibration\My_CCW'
    settings={'dir': workingDirectory,
                'x': 'Ch2_My', 'y': 'Moment',
                'xUnit': '[volts]', 'yUnit': '[Newton Meter]',
                'xUncert': 0.0, 'yUncert': np.sqrt((armLength*weightUncert)**2 + (np.cos(momentArmTilt)*maxWeight*armLengthUncert)**2),
                'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botMCCW=Calibrator(**settings)
    botMCCW.read_data()
    discrete_fit_data_MCCW=decorator(botMCCW, botMCCW.discrete_fit_data, h=lambda y: y/1000*9.81*np.cos(momentArmTilt)*armLength)
    resultMCCW=discrete_fit_data_MCCW(customEquation)
    MomentCCWData=[
            ['Parameter', 'Parameter_Std', 'R2', 'Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultMCCW['param'][0], resultMCCW['paramStd'][0], resultMCCW['R2'], max(resultMCCW['scatterUncert']), max(resultMCCW['meanUncert'])],
            [resultMCCW['param'][1], resultMCCW['paramStd'][1], '', '', '']
            ]

    # NOTE: CW Moment Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\WS16_VWT-V10-01_calibration\My_CW'
    settings={'dir': workingDirectory,
                'x': 'Ch2_My', 'y': 'Moment',
                'xUnit': '[volts]', 'yUnit': '[Newton Meter]',
                'xUncert': 0.0, 'yUncert': np.sqrt((armLength*weightUncert)**2 + (np.cos(momentArmTilt)*maxWeight*armLengthUncert)**2),
                'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botMCW=Calibrator(**settings)
    botMCW.read_data()
    discrete_fit_data_MCW=decorator(botMCW, botMCW.discrete_fit_data, h=lambda y: -y/1000*9.81*np.cos(momentArmTilt)*armLength)
    resultMCW=discrete_fit_data_MCW(customEquation)
    MomentCWData=[
            ['Parameter', 'Parameter_Std', 'R2', 'Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultMCW['param'][0], resultMCW['paramStd'][0], resultMCW['R2'], max(resultMCW['scatterUncert']), max(resultMCW['meanUncert'])],
            [resultMCW['param'][1], resultMCW['paramStd'][1], '', '', '']
            ]

    # NOTE: Total Moment Calibration
    workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\WS16_VWT-V10-01_calibration\My'
    settings={'dir': workingDirectory,
                'x': 'Ch2_My', 'y': 'Moment',
                'xUnit': '[volts]', 'yUnit': '[Newton Meter]',
                'xUncert': 0.0, 'yUncert': np.sqrt((armLength*weightUncert)**2 + (np.cos(momentArmTilt)*maxWeight*armLengthUncert)**2),
                'auto': True,
                'eqn': 'p1'}
    customEquation=None
    botM=Calibrator(**settings)
    botM.read_data()
    discrete_fit_data_M=decorator(botM, botM.discrete_fit_data, h=lambda y: y/1000*9.81*np.cos(momentArmTilt)*armLength)
    resultMoment=discrete_fit_data_M(customEquation)
    MomentData=[
            ['Parameter', 'Parameter_Std', 'R2', 'Max95%Scatter_Uncertainty', 'Max95%Mean_Uncertainty'],
            [resultMoment['param'][0], resultMoment['paramStd'][0], resultMoment['R2'], max(resultMoment['scatterUncert']), max(resultMoment['meanUncert'])],
            [resultMoment['param'][1], resultMoment['paramStd'][1], '', '', '']
            ]

    os.chdir(savingDirectory)
    calibrationResult=collections.OrderedDict()
    calibrationResult.update({'SG': SGData, 'AOA': AOAData, 'Drag': DragData,
                        'MomentCCW': MomentCCWData, 'MomentCW': MomentCWData,
                        'MomentTotal': MomentData})
    pyexcel.isave_book_as(bookdict=calibrationResult, dest_file_name='results.xlsx')


# velocity=11.25
# velocityUncertanties=0.05
# density=1.184
# densityUncertanties=0.2/100*density
# chord=0.14
# chordUncertanties=0
# span=0.197
# spanUncertanties=0
# airfoilMass=149.8/1000
#
# fSG=resultSG['f']
# fAOA=resultAOA['f']
# fDrag=resultDrag['f']
# fMoment=resultMoment['f']

# NOTE: CL vs AOA
# CL=lambda y, v, rho, c, s: 2*y/(rho*c*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.2.2018\Re 100 000\target_AOA'
# settings={'dir': workingDirectory,
#             'x': 'Ch4_a_control', 'y': 'Ch1_Fz',
#             'xUnit': '[degree]', 'yUnit': '[.]',
#             'xUncert': levelUncert, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Re100000_targetAOA1.csv',
#             'eqn': 'p1'}
# customEquation=None
# botCL=Calibrator(**settings)
# inputs={'xStd': max(resultAOA['scatterUncert']), 'yStd': max(resultSG['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':1}
# botCL.read_data()
# analyzerCL=decorator(botCL, botCL.analyzer,
#         g=lambda x:fAOA(x,*resultAOA['param']), h=lambda y:fSG(y,*resultSG['param'])+9.81*airfoilMass)
# resultCL=analyzerCL(fx=lambda x:x, fy=CL, inputs=inputs)
# os.chdir(savingDirectory)
# CLResults=collections.OrderedDict()
# [CLResults.update({key: [list(resultCL[key].keys()), *np.array(list(resultCL[key].values())).transpose().tolist()]}) for key in resultCL.keys()]
# pyexcel.isave_book_as(bookdict=CLResults, dest_file_name='AOA_results_CL_Re_1E5.xlsx')

# NOTE: CD vs AOA
# CD=lambda y, v, rho, c, s: 2*y/(rho*c*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.2.2018\Re 100 000\target_AOA'
# settings={'dir': workingDirectory,
#             'x': 'Ch4_a_control', 'y': 'Ch0_Fx',
#             'xUnit': '[degree]', 'yUnit': '[.]',
#             'xUncert': levelUncert, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Re100000_targetAOA1.csv',
#             'eqn': 'p1'}
# customEquation=None
# botCD=Calibrator(**settings)
# inputs={'xStd': max(resultAOA['scatterUncert']), 'yStd': max(resultDrag['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':1}
# botCD.read_data()
# analyzerDrag=decorator(botCD, botCD.analyzer,
#         g=lambda x:fAOA(x,*resultAOA['param']), h=lambda y:fDrag(y,*resultDrag['param']))
# resultCD=analyzerDrag(fx=lambda x:x, fy=CD, inputs=inputs)
# os.chdir(savingDirectory)
# CDResults=collections.OrderedDict()
# [CDResults.update({key: [list(resultCD[key].keys()), *np.array(list(resultCD[key].values())).transpose().tolist()]}) for key in resultCD.keys()]
# pyexcel.isave_book_as(bookdict=CDResults, dest_file_name='AOA_results_CD_Re_1E5.xlsx')

# NOTE: CM vs AOA
# CM=lambda y, v, rho, c, s: 2*y/(rho*c**2*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.2.2018\Re 100 000\target_AOA'
# settings={'dir': workingDirectory,
#             'x': 'Ch4_a_control', 'y': 'Ch2_My',
#             'xUnit': '[degree]', 'yUnit': '[.]',
#             'xUncert': levelUncert, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Re100000_targetAOA1.csv',
#             'eqn': 'p1'}
# customEquation=None
# botCM=Calibrator(**settings)
# inputs={'xStd': max(resultAOA['scatterUncert']), 'yStd': max(resultMoment['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':1}
# botCM.read_data()
# analyzerMoment=decorator(botCM, botCM.analyzer,
#         g=lambda x: fAOA(x, *resultAOA['param']), h=lambda y:fMoment(y, *resultMoment['param']))
# resultCM=analyzerMoment(fx=lambda x:x, fy=CM, inputs=inputs)
# os.chdir(savingDirectory)
# CMResults=collections.OrderedDict()
# [CMResults.update({key: [list(resultCM[key].keys()), *np.array(list(resultCM[key].values())).transpose().tolist()]}) for key in resultCM.keys()]
# pyexcel.isave_book_as(bookdict=CMResults, dest_file_name='AOA_results_CM_Re_1E5.xlsx')


# Fog Machine
# velocity=1.5
# velocityUncertanties=0.05
# density=1.188
# densityUncertanties=0.2/100*density
# chord=0.14
# chordUncertanties=0
# span=0.197
# spanUncertanties=0
# airfoilMass=149.8/1000
#
# fSG=resultSG['f']
# fAOA=resultAOA['f']
# fDrag=resultDrag['f']
# fMoment=resultMoment['f']


# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.9.2018'
# settings={'dir': workingDirectory,
#             'x': 'X_Value', 'y': 'Ch4_a_control',
#             'xUnit': '[Second]', 'yUnit': '[degree]',
#             'xUncert': 0.0, 'yUncert': levelUncert,
#             'auto': True, 'target_file': 'Fog_Richard_Liu_Sweep.csv',
#             'eqn': 'p1'}
# customEquation=None
# botFogAOA=Calibrator(**settings)
# inputs={'xStd': 0, 'yStd': max(resultAOA['scatterUncert']),
#         'fxParam': [], 'fyParam': [],
#         'fxUncert': [], 'fyUncert': [],
#         'threshold':2000}
# botFogAOA.read_data()
# analyzerFogAOA=decorator(botFogAOA, botFogAOA.analyzer,
#                         h=lambda y:fAOA(y,*resultAOA['param']))
# resultFogAOA=analyzerFogAOA(fx=lambda x:x, fy=lambda y:y, inputs=inputs)
# os.chdir(savingDirectory)
# FogAOAResults=collections.OrderedDict()
# [FogAOAResults.update({key: [list(resultFogAOA[key].keys()), *np.array(list(resultFogAOA[key].values())).transpose().tolist()]}) for key in resultFogAOA.keys()]
# pyexcel.isave_book_as(bookdict=FogAOAResults, dest_file_name='Fog_Time_AOA_Re_1E4.xlsx')



# CL=lambda y, v, rho, c, s: 2*y/(rho*c*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.9.2018'
# settings={'dir': workingDirectory,
#             'x': 'X_Value', 'y': 'Ch1_Fz',
#             'xUnit': '[Second]', 'yUnit': '[.]',
#             'xUncert': 0.0, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Fog_Richard_Liu_Sweep.csv',
#             'eqn': 'p1'}
# customEquation=None
# botFogCL=Calibrator(**settings)
# inputs={'xStd': 0, 'yStd': max(resultSG['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':2000}
# botFogCL.read_data()
# analyzerFogCL=decorator(botFogCL, botFogCL.analyzer,
#                         h=lambda y:fSG(y,*resultSG['param'])+9.81*airfoilMass)
# resultFogCL=analyzerFogCL(fx=lambda x:x, fy=CL, inputs=inputs)
# os.chdir(savingDirectory)
# FogCLResults=collections.OrderedDict()
# [FogCLResults.update({key: [list(resultFogCL[key].keys()), *np.array(list(resultFogCL[key].values())).transpose().tolist()]}) for key in resultFogCL.keys()]
# pyexcel.isave_book_as(bookdict=FogCLResults, dest_file_name='Fog_results_CL_Re_1E4.xlsx')


# CD=lambda y, v, rho, c, s: 2*y/(rho*c*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.9.2018'
# settings={'dir': workingDirectory,
#             'x': 'X_Value', 'y': 'Ch0_Fx',
#             'xUnit': '[Second]', 'yUnit': '[.]',
#             'xUncert': 0.0, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Fog_Richard_Liu_Sweep.csv',
#             'eqn': 'p1'}
# customEquation=None
# botFogCD=Calibrator(**settings)
# inputs={'xStd': 0, 'yStd': max(resultDrag['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':2000}
# botFogCD.read_data()
# analyzerFogCD=decorator(botFogCD, botFogCD.analyzer,
#                                     h=lambda y:fDrag(y,*resultDrag['param']))
# resultFogCD=analyzerFogCD(fx=lambda x:x, fy=CD, inputs=inputs)
# os.chdir(savingDirectory)
# FogCDResults=collections.OrderedDict()
# [FogCDResults.update({key: [list(resultFogCD[key].keys()), *np.array(list(resultFogCD[key].values())).transpose().tolist()]}) for key in resultFogCD.keys()]
# pyexcel.isave_book_as(bookdict=FogCDResults, dest_file_name='Fog_results_CD_Re_1E4.xlsx')


# CM=lambda y, v, rho, c, s: 2*y/(rho*c**2*s*v**2)
# workingDirectory=r'C:\Users\Richard\Google Drive\Berkeley\Berkeley Classes\Fall 2018\ME 103\Lab 4-5\ME103\11.9.2018'
# settings={'dir': workingDirectory,
#             'x': 'X_Value', 'y': 'Ch2_My',
#             'xUnit': '[Second]', 'yUnit': '[.]',
#             'xUncert': 0.0, 'yUncert': 0.0,
#             'auto': True, 'target_file': 'Fog_Richard_Liu_Sweep.csv',
#             'eqn': 'p1'}
# customEquation=None
# botFogCM=Calibrator(**settings)
# inputs={'xStd': 0, 'yStd': max(resultMoment['scatterUncert']),
#         'fxParam': [], 'fyParam': [velocity, density, chord, span],
#         'fxUncert': [], 'fyUncert': [velocityUncertanties, densityUncertanties, chordUncertanties, spanUncertanties],
#         'threshold':2000}
# botFogCM.read_data()
# analyzerFogCM=decorator(botFogCM, botFogCM.analyzer,
#                                 h=lambda y:fMoment(y, *resultMoment['param']))
# resultFogCM=analyzerFogCM(fx=lambda x:x, fy=CM, inputs=inputs)
# os.chdir(savingDirectory)
# FogCMResults=collections.OrderedDict()
# [FogCMResults.update({key: [list(resultFogCM[key].keys()), *np.array(list(resultFogCM[key].values())).transpose().tolist()]}) for key in resultFogCM.keys()]
# pyexcel.isave_book_as(bookdict=FogCMResults, dest_file_name='Fog_results_CM_Re_1E4.xlsx')





























plt.show()
plt.close()
