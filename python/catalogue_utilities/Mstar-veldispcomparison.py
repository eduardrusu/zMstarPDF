# Cone search
import numpy as np
import pylab as plt

# section 3, Jabran Zahid et al. 2016
logsigmab = 2.072
logMb = 10.26
alpha1 = 0.403
alpha2 = 0.29

def log_sigma_Zahid(logMstar):
    log_sigma_Zahid = np.copy(logMstar)
    log_sigma_Zahid[logMstar <= logMb] = logsigmab + alpha1*(logMstar[logMstar <= logMb] - logMb)
    log_sigma_Zahid[logMstar > logMb] = logsigmab + alpha2*(logMstar[logMstar > logMb] - logMb)
    return log_sigma_Zahid

logMstar = np.linspace(9,12,100)

# bottom of left column, page 5, Mason et al. 2015
p = 0.24
q = 2.34
log_sigma_Mason = p * logMstar - 11 * p + q


plt.plot(logMstar,log_sigma_Mason,label="Mason 2015")
plt.plot(logMstar,log_sigma_Zahid(logMstar),label="Zahid 2016")
plt.xlabel('log(Mstar)')
plt.ylabel('log(sigma)')
plt.legend()
plt.show()
