import L_R_zams

Z = input('Enter Z: ')
dzeta = log(Z/0.02)


data1 = np.loadtxt("coeffMS.txt", delimiter='\t', dtype=np.float)
a = []
for i in range(34):
	a.append(data1[i][0] + data1[i][1]*dzeta + data1[i][2]*dzeta**2 + data1[i][3]*dzeta**3 + data1[i][4]*dzeta**4)

M_hook = 1.0185 + 0.16015*dzeta + 0.0892*dzeta**2
M_HeF = 1.995 + 0.25*dzeta + 0.087*dzeta**2
M_FGB = (13.048*(Z/0.002)**0.06)/(1 + 0.0012*(0.002/Z)**1.27)
#???M_CE = (1.0/3)*M_E

x = max(0.95,min(0.95 - 0.03*(dzeta + 0.30103), 0.99))
mu = max(0.5, 1.0 - 0.01*max(a[5]/M**a[6] , a[7] + a[8]/M**a[9]))
t_BGB = (a[0] + a[1]*M**4 + a[2]*M**5.5 + M**7)/(a[3]*M**2 + a[4]*M**7)
t_hook = mu*t_BGB
t_MS = max(t_hook, x*t_BGB)

#L_ZAMS
#R_ZAMS
a26 = 5.5
c1 = -0.08672073 
L_TMS = (a[10]*M**3 + a[11]*M**4 + a[12]*M**9)/(a[13] + a[14]*M**5 + M**7.2)
R_TMS = (c1*M**3 + a[16]*M**a26 +a[17]*M**(a26+1.5))/(a[18] + M**5)

def MS(M, t, t_MS, L_ZAMS, R_ZAMS):
	t = t + 0.1
	eps = 0.01
	tau = t/t_MS
	tau1 = min(1.0, t/t_hook)
	tau2 = max(0.0, min(1.0, (t - (1.0 - eps)*t_hook)/(eps*t_hook)))
	
	B2 = a76 + a77*(1.0 - a78)**a79
	C = a80

	if Z >= 0.0009:
		eta = 10
	else:
		eta = 20

	a71 = 3.5
	if 2.0 <= M <= 16.0:
		betaRi = (a[31]*M**3.5)/(a[32] + M**a71)
	elif 16.0 < M:
		betaRi = (a[31]*16.0**3.5)/(a[32] + 16**a71) + a[33]*(M - 16.0) 
	betaR = betaRi - 1.0

	a37 = 0.6
	a41 = 3.6
	a48 = 1.6
	a57 = 1.5
	a56 = 0.96
	a60 = 1.8
	a61 = 2.3
	a67 = 5
	gamma = 0
	C1 = (a[28]*M**a60)/(a[29] + M**a61)
	B1 = max(0.0, a[34] - a[35]*a57**a56)
	aL = (a[25] + a[26]*M**a48)/(M**0.4 + a[27]*M**1.9)
	aR = C1 + a[30]*(M - a67)
	betaL = max(0.0, B1 - 10.0*(M - a57)*B1)
	deltaL = min(a[19]/(M**a[20]), a[21]/(M**a37))
	deltaR = (a[22] + a[23]*M**3.5)/(a[24]*M**3 + M**a41) - 1.0

	#L_b = exp(deltaL*(tau1**2 - tau2**2))
	#L_MS = L_a/L_b

	log(L_MS/L_ZAMS) = aL*tau + betaL*tau**eta + (log(L_TMS/L_ZAMS) - aL - betaL)*tau**2 - deltaL(tau1**2 - tau2**2)
	log(R_MS/R_ZAMS) = aR*tau _ betaR*tau**10 + gamma*tau**40 + (log(R_TMS/R_ZAMS) - aR - betaR - gamma)*tau**3 - deltaR(tau1**3 - tau2**3)