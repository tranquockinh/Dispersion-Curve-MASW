from Paras import *
def MASW_Despersion_image(u,N,fs,cT_min,cT_max,delta_cT):
    
    ## Converting measuring frequency from Hz to rad/sec
    omega_fs = 2 * np.pi * fs
    ## Number of samples in each trace
    Lu = len(u[:,0])
    ## Empty matrices with Lu lines and n columns
    U = np.zeros((Lu,N),dtype='complex')
    P = np.zeros((Lu,N),dtype='complex')
    Unorm = np.zeros((Lu,N),dtype='complex')
    
    ## Waveform decomposition (Fourier transform)
    for j in range(N):    
        U[:,j] = np.fft.fft(u[:,j])
    
    LU = len(U[:,0])
    ## Normalize U in offset and frequency domain
    ## Compute the phase spectrum of UnboundLocalError
    for j in range(N):
        for k in range(int(LU)):
            Unorm[k,j] = U[k,j] / np.abs(U[k,j])
        P[:,j] = np.exp(1j * (-np.angle(U[:,j])))
    ## Frequency range for U
    omega = (1/LU) * np.arange(0,LU+1,1) * omega_fs
    ## compute Slant-stack (summed) amplitude corresponding to
    ## each set of [omega, cT, A(omega,cT)]
    cT = np.arange(cT_min, cT_max + delta_cT, delta_cT)
    LcT = len(cT)
    
    ## Empty matrices with Lu lines and n columns
    c = np.zeros((LU, LcT))
    f = np.zeros((LU,LcT))
    A = np.zeros((LU,LcT))
    
    for j in range(int(LU)):
        for k in range(int(LcT)):
        ## Frequency (in [Hz]) corresponding to angular frequency omega
            f[j,k] = omega[j] / (2 * np.pi)
            ## Testing phase velocity [m/s]
            c[j,k] = cT[k]
            ## Determining the amount of phase shifts required to counterbalance
            ## the time delay corresponding to specific offsets for a given set 
            ## of omega and cT
            delta = omega[j] / cT[k]
            ## Applying the phase shifts (delta) to distinct traces of the 
            ## transformed shot gather
            ## Obtaining the (normalized) slant-stack (summed) amplitude
            ## corresponding to each set of omega and cT
            temp = 0
            for l in range(N):
                temp = temp + np.exp(-1j * delta * x[l]) * P[j,l]
            
            ## Compute absolute value and normalize with respect to number of 
            ## receiver
            A[j,k] = np.abs(temp) / N
    
    return f,c,A