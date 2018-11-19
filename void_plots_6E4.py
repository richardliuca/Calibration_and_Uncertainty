import matplotlib.pyplot as plt
import numpy as np
# NOTE: Re:6E4
# NOTE: Experimental
AOA=[-25.10,-20.72,-15.82,-11.43,-6.20,-1.72,-1.40,1.59,4.26,6.49,8.20,10.10,
11.78,13.58,15.07,16.80,18.36,20.12,23.69]
AOA_E=[0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,
0.24,0.24,0.24,0.24,0.24]

CL=[-0.56,-0.49,-0.45,-0.46,-0.36,-0.16,-0.19,-0.01,0.20,0.31,0.40,0.50,0.56,
0.61,0.63,0.55,0.48,0.49,0.49]
CL_E=[0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,
0.12,0.12,0.12,0.12,0.12]

CD=[0.33,0.25,0.17,0.12,0.03,0.01,0.02,0.01,0.02,0.02,0.03,0.04,0.05,0.06,0.08,
0.16,0.19,0.21,0.24]
CD_E=[0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,
0.08,0.08,0.08,0.08,0.08]

CM=[-0.35,-0.34,-0.35,-0.34,-0.29,-0.30,-0.31,-0.31,-0.30,-0.31,-0.31,-0.33,
-0.33,-0.33,-0.33,-0.27,-0.25,-0.24,-0.23]
CM_E=[0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,
0.15,0.15,0.15,0.15,0.15]
# NOTE: CFD
A_CFD=[0,4,8,12]
CL_CFD=[0.261,0.734,1.139,1.144]
CD_CFD=[0.01532,0.01841,0.02794,0.10236]
CM_CFD=[-0.05,-0.055,-0.059,-0.027]


fig, axs=plt.subplots(3,1,sharex=True)
axs[0].plot(AOA, CL, label='Experimental')
axs[0].errorbar(AOA, CL, xerr=AOA_E, yerr=CL_E, fmt='.', label='Scatter Uncertanty' )
axs[0].plot(A_CFD, CL_CFD, label='CFD')
axs[0].set_ylabel('CL',size=26)
axs[0].grid(b=True, which='both', axis='both')
axs[0].legend(fontsize=16)
axs[0].tick_params(labelsize=18)

axs[1].plot(AOA, CD, label='Experimental')
axs[1].errorbar(AOA, CD, xerr=AOA_E, yerr=CD_E, fmt='.', label='Scatter Uncertanty' )
axs[1].plot(A_CFD, CD_CFD, label='CFD')
axs[1].set_ylabel('CD',size=26)
axs[1].grid(b=True, which='both', axis='both')
axs[1].legend(fontsize=16)
axs[1].tick_params(labelsize=18)

axs[2].plot(AOA, CM, label='Experimental')
axs[2].errorbar(AOA, CM, xerr=AOA_E, yerr=CM_E, fmt='.', label='Scatter Uncertanty' )
axs[2].plot(A_CFD, CM_CFD, label='CFD')
axs[2].grid(b=True, which='both', axis='both')
axs[2].legend(fontsize=16)
axs[2].set_ylabel('CM',size=26)
axs[2].set_xlabel('Angle of Attack [degree]',size=26)
axs[2].tick_params(labelsize=18)


plt.draw()
plt.show()
