import matplotlib.pyplot as plt
import numpy as np
# NOTE: Re:1E5
# NOTE: Experimental
AOA=[-25.00,-21.17,-11.90,-10.14,-6.98,-1.78,1.43,3.72,6.15,7.96,10.13,11.36,
13.13,14.05,15.69,17.31,18.28,19.56,23.68]
CL=[-0.43,-0.39,-0.36,-0.35,-0.25,-0.04,0.15,0.28,0.37,0.45,0.53,0.56,
0.58,0.55,0.49,0.44,0.44,0.45,0.46]
CD=[0.25,0.19,0.09,0.03,0.02,0.01,0.01,0.02,0.03,0.03,0.04,0.05,0.07,0.09,0.15,
0.17,0.19,0.21,0.22]
CM=[-0.14,-0.13,-0.12,-0.07,-0.08,-0.10,-0.10,-0.11,-0.12,-0.12,-0.13,-0.13,
-0.13,-0.12,-0.07,-0.06,-0.06,-0.06,-0.06]

# NOTE: CFD
A_CFD=[0,4,8,12]
CL_CFD=[0.261,0.733,1.139,1.144]
CD_CFD=[0.01197,0.01483,0.02418,0.09473]
CM_CFD=[-0.051,-0.055,-0.059,-0.029]


fig, axs=plt.subplots(3,1,sharex=True)
axs[0].plot(AOA, CL, label='Experimental')
axs[0].plot(A_CFD, CL_CFD, label='CFD')
axs[0].set_ylabel('CL',size=26)
axs[0].grid(b=True, which='both', axis='both')
axs[0].legend(fontsize=16)
axs[0].tick_params(labelsize=12)

axs[1].plot(AOA, CD, label='Experimental')
axs[1].plot(A_CFD, CD_CFD, label='CFD')
axs[1].set_ylabel('CD',size=26)
axs[1].grid(b=True, which='both', axis='both')
axs[1].legend(fontsize=16)
axs[1].tick_params(labelsize=18)

axs[2].plot(AOA, CM, label='Experimental')
axs[2].plot(A_CFD, CM_CFD, label='CFD')
axs[2].grid(b=True, which='both', axis='both')
axs[2].legend(fontsize=16)
axs[2].set_ylabel('CM',size=26)
axs[2].set_xlabel('Angle of Attack [degree]',size=26)
axs[2].tick_params(labelsize=18)


plt.draw()
plt.show()
