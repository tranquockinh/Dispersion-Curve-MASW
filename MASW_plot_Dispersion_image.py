from Paras import *
def MASW_plot_Dispersion_image(f,c,A):
    
    resolution = 100
    fmin = 1 ## Hz
    fmax = 50 ## Hz
    ## Limit of frequency axis
    RemoveMin = np.abs((f[:,0] - fmin))
    RemoveMax = np.abs((f[:,0] - fmax))
    IdxMin = np.array(np.where(RemoveMin==np.min(RemoveMin))).flatten()[0]
    IdxMax = np.array(np.where(RemoveMax==np.min(RemoveMax))).flatten()[0]
    valMin, valMax = RemoveMin[IdxMin], RemoveMax[IdxMax]
    Aplot = A[IdxMin:IdxMax+1,:]
    fplot = f[IdxMin:IdxMax+1,:]
    cplot = c[IdxMin:IdxMax+1,:]
    
    #fig, ax = plt.subplots(dpi=100)
    #plt.contourf(fplot, cplot, Aplot, resolution, cmap=plt.cm.jet)    
    #ax.set_xticks(np.arange(0,10+1E-5,fmax+0.01))
    #ax.set_xlim(fmin,fmax)
    #ax.set_xlabel('Frequency (Hz)', fontsize = FigFontSize)
    #ax.set_ylabel('Phase velocity [m/s]', fontsize = FigFontSize)
    
    #cbar = plt.colorbar(orientation='horizontal',location='top')
    #cbar.set_label('Normalized amplitude', fontsize = FigFontSize)
    #plt.show()
    
    return valMin, valMax, IdxMin, IdxMax, fplot, cplot, Aplot