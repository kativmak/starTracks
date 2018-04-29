import numpy as np
from math import log10, log
import matplotlib.pyplot as plt
 
data1 = np.loadtxt("table1.txt", delimiter='\t', dtype=np.float)
data2 = np.loadtxt("table2.txt", delimiter='\t', dtype=np.float)
lettersL = []
lettersR = []
Z_sun = 0.02
Z = 0.001
a = 0.1
ratio = log10(Z/Z_sun)

for i in range(7):
	lettersL.append(data1[i][0] + data1[i][1]*ratio + data1[i][2]*ratio**2 + data1[i][3]*ratio**3 + data1[i][4]*ratio**4)
	
for i in range(9):
	lettersR.append(data2[i][0] + data2[i][1]*ratio + data2[i][2]*ratio**2 + data2[i][3]*ratio**3 + data2[i][4]*ratio**4)

M = np.linspace(0.1, 100, 601)
L = lambda M: ((lettersL[0]*M**5.5 + lettersL[1]*M**11)/
	(lettersL[2] + M**3 + lettersL[3]*M**5 + lettersL[4]*M**7 + lettersL[5]*M**8 + lettersL[6] *M**9.5))
R = lambda M: ((lettersR[0]*M**2.5 + lettersR[1]*M**6.5 + lettersR[2]*M**11 + lettersR[3]*M**19 + lettersR[4]*M**19.5)/
	(lettersR[5] + lettersR[6]*M**2 + lettersR[7]*M**8.5 + M**18.5 + lettersR[8]*M**19.5))

l = []
r = []
Mlg = []
for index, item in enumerate(M):
	l.append(log10(L(item)))
	r.append(log10(R(item)))
	Mlg.append(log(item))

plt.plot(Mlg, l , color = 'k', linestyle ='--')
plt.xlabel(r'$log(M/M_{sun})$')
plt.ylabel(r'$log(L/L_{sun})$')
plt.grid(True)
plt.show()

#save(name='pic_2_1', fmt='pdf')