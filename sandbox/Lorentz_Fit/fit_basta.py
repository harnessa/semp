import numpy as np
import matplotlib.pyplot as plt;plt.ion()

# einf =        #infinite epsilon for silicon
# eps0 =        #permitivitty of freespace
# e =           #charge of electron
# me =          #mass of electron
# eV_um_scale = um_scale/1.23984193

def A_func(angfreq, gamma, beta, mu, omega):
    return np.exp(1j*beta) * (omega - angfreq - 1j*gamma)**mu

def B_func(angfreq, gamma, beta, mu, omega):
    return np.exp(-1j*beta) * (omega + angfreq + 1j*gamma)**mu

def C_func(angfreq, gamma, beta, mu, omega):
    return 2*np.real(np.exp(-1j*beta) * (omega + 1j*gamma)**mu)

def D_func(angfreq, gamma, beta, mu, omega):
    return 2*1j*mu*angfreq * np.imag(np.exp(-1j*beta) * (omega + 1j*gamma)**(mu-1))

def ABCD_func(angfreq, gbmo):
    return A_func(angfreq, *gbmo) + B_func(angfreq, *gbmo) - \
        C_func(angfreq, *gbmo) - D_func(angfreq, *gbmo)

def E_func(angfreq, Neh, tD, mopt):
    return Neh * e**2 / (eps0 * mopt * me * angfreq**2) / (1 + 1j/(tD*angfreq))


def epsilon(angfreq, params):
    ans = einf - E_func(angfreq)
    for pms in params:
        ans += pms[0]/angfreq**2 * ABCD_func(angfreq, pms[1:])
    return ans

#Nd = 1e14
cgbme = np.array([[69.56, 239.3, 92], [0.09775, 0.361, 0.2316], [0.3582, 0.307, 0.004174],
    [0.6976, 0.4398, 1.141], [3.368, 3.654, 4.287]]).T


#Convert critical energy to angfreq
# cgbme[:,-1] =

breakpoint()
