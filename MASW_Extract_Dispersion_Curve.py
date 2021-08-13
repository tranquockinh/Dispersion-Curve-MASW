from Paras import *
from MASW_plot_Dispersion_image import MASW_plot_Dispersion_image

def MASW_Extract_Dispersion_Curve(f,c,A):
    resolution = 100
    fmin = 1 ## Hz
    fmax = 50 ## Hz
    f_receivers = 4.5 ## Hz
    up_low_boundary = 'yes' ## The region, points inside supposed belonged to Fundamental DC
    p = 95 ## Percentage value for determination of upper/lower
    resolution = 100
    
    valMin, valMax, IdxMin, IdxMax, fplot, cplot, Aplot  = MASW_plot_Dispersion_image(f,c,A)
    fig, ax = plt.subplots(dpi=100)
    plt.contourf(fplot, cplot, Aplot, resolution, cmap=plt.cm.jet)    
    ax.set_xticks(np.arange(0,10+1E-5,fmax+0.01))
    ax.set_xlim(fmin,fmax)
    ax.set_xlabel('Frequency (Hz)', fontsize = FigFontSize)
    ax.set_ylabel('Phase velocity [m/s]', fontsize = FigFontSize)
    cbar = plt.colorbar(orientation='horizontal',location='top')
    cbar.set_label('Normalized amplitude', fontsize = FigFontSize)
    
    Aabsnorm2 = np.zeros((Aplot.shape))
    
    for i in range(len(Aplot[:,0])):
        Aabsnorm2[i,:] = np.divide(Aplot[i,:] , np.max(Aplot[i,:]))
    
    loc = np.argwhere(Aabsnorm2.T==1)
    ## Reurn the first column is the search results on y-axis
    ## Rturn the second column is the search results on x-axis
    f_loc, c_loc = loc[:,1], loc[:,0]
    
    Amax_fvec = np.zeros((len(f_loc)))
    Amax_cvec = np.zeros((len(c_loc)))
        
    for i in range(len(f_loc)):
        Amax_fvec[i] = fplot[f_loc[i], 0]
        Amax_cvec[i] = cplot[0, c_loc[i]]
        
    ## Extract those values of freq. with max amplitude > freq. of receivers
    ii = np.where(Amax_fvec > f_receivers)
    Amax_fvec = Amax_fvec[ii]
    Amax_cvec = Amax_cvec[ii]
    
    ## Sort freq. vectore of max amplitude in order
    jj = np.argsort(Amax_fvec) ## return list of sorted indices
    Amax_fvec_sort = np.zeros((len(jj)))
    Amax_cvec_sort = np.zeros((len(jj)))  
    '''
    for i, j in enumerate(jj):
        Amax_fvec_sort[i] = Amax_fvec[j]
        Amax_cvec_sort[i] = Amax_cvec[j]
    '''
    Amax_fvec_sort = Amax_fvec[jj] 
    Amax_cvec_sort = Amax_cvec[jj]
        
    if up_low_boundary == 'yes':
        
        loc_p = np.argwhere(Aabsnorm2.T > p/100)
        f_loc_p, c_loc_p = loc_p[:,1], loc_p[:,0]
        
        Amax_fvec_p = np.zeros((f_loc_p.shape))
        Amax_cvec_p = np.zeros((c_loc_p.shape))
        
        for i in range(len(f_loc_p)):                
            Amax_fvec_p[i] = fplot[f_loc_p[i], 0]
            Amax_cvec_p[i] = cplot[0,c_loc_p[i]]
        
        ## Search higher reveiver frequenct components
        kk = np.where(Amax_fvec_p > f_receivers)
        Amax_fvec_p = Amax_fvec_p[kk]
        Amax_cvec_p = Amax_cvec_p[kk]
        
        pp = np.argsort(Amax_fvec_p)
        Amax_fvec_sort_p = Amax_fvec_p[pp]
        Amax_cvec_sort_p = Amax_cvec_p[pp]
        
        ## Compute up-low-boundary points
        Amax_fvec_sort_p_cell = np.empty((len(np.unique(Amax_fvec_sort_p)), 1),dtype=object)
        Amax_cvec_sort_p_cell = np.empty((len(np.unique(Amax_fvec_sort_p)), 1),dtype=object)
        f_curve0_up_temp = np.zeros((len(np.unique(Amax_fvec_sort_p)),1))
        c_curve0_up_temp = np.zeros((len(np.unique(Amax_fvec_sort_p)),1))
        f_curve0_low_temp = np.zeros((len(np.unique(Amax_fvec_sort_p)),1))
        c_curve0_low_temp = np.zeros((len(np.unique(Amax_fvec_sort_p)),1))
        
        U = np.unique(Amax_fvec_sort_p, return_index=False) ## return index's sorted order
        Unique_Amax_fvec_sort_p = np.empty((len(U),1),dtype=object)
        
        for i in range(len(U)):
                        
            ## Store each ndarray in a cell (dtype = object)
            Amax_fvec_sort_p_cell[i][0] = Amax_fvec_sort_p[np.where(Amax_fvec_sort_p==U[i])]
            Amax_cvec_sort_p_cell[i][0] = Amax_cvec_sort_p[np.where(Amax_fvec_sort_p==U[i])]
        
            f_curve0_up_temp[i] = np.max(Amax_fvec_sort_p_cell[i][0])
            c_curve0_up_temp[i] = np.max(Amax_cvec_sort_p_cell[i][0])
            f_curve0_low_temp[i] = np.min(Amax_fvec_sort_p_cell[i][0])
            c_curve0_low_temp[i] = np.min(Amax_cvec_sort_p_cell[i][0])
    
    # Plotting
    fig, ax = plt.subplots(dpi=100)
    plt.contourf(fplot, cplot, Aplot, resolution, cmap=plt.cm.jet)    
    #ax.set_xticks(np.arange(0,10+1E-5,fmax+0.01))
    #ax.set_xlim(fmin,fmax)
    ax.set_xlabel('Frequency (Hz)', fontsize = FigFontSize)
    ax.set_ylabel('Phase velocity [m/s]', fontsize = FigFontSize)
    
    cbar = plt.colorbar(orientation='horizontal',location='top')
    cbar.set_label('Normalized amplitude', fontsize = FigFontSize)
    
    plt.plot(Amax_fvec_sort, Amax_cvec_sort, 'o',\
             markersize=4,markerfacecolor='k',markeredgecolor='k')   
    
    plt.plot(Amax_fvec_sort_p,Amax_cvec_sort_p,'o',\
            markersize=1,markerfacecolor='k',markeredgecolor='k')
    ## Plot up-low boundary
    
    plt.plot(f_curve0_up_temp,c_curve0_up_temp, 'o',\
             markersize=4,markerfacecolor='k',markeredgecolor='k')
    plt.plot(f_curve0_low_temp, c_curve0_low_temp,'o',\
             markersize=4,markerfacecolor='k',markeredgecolor='k')
    
    labels = np.arange(1,len(Amax_fvec_sort)+1,1)
    
    for i, addtext in enumerate(labels):
        ax.annotate(addtext, (Amax_fvec_sort[i], Amax_cvec_sort[i]),\
        horizontalalignment='right', verticalalignment='bottom')
        
    ##plt.show()
    
    lower_freq = int((input('Input lower frequency value: ')))
    higher_freq = int((input('Input higher frequency value: ')))
    nP0 = np.arange(lower_freq, higher_freq+1,1)
    f_curve0 = Amax_fvec_sort[nP0]
    c_curve0 = Amax_cvec_sort[nP0]
    
    if up_low_boundary == 'yes':
        f_curve0_up = f_curve0_up_temp[nP0]
        c_curve0_up = c_curve0_up_temp[nP0]
        f_curve0_low = f_curve0_low_temp[nP0]
        c_curve0_low = c_curve0_low_temp[nP0]
    
    lambda_curve0 = np.divide(c_curve0, f_curve0)
    
    if up_low_boundary == 'yes':
        lambda_curve0_up = np.divide(c_curve0_up, f_curve0_up)
        lambda_curve0_low = np.divide(c_curve0_low, f_curve0_low)
    return  f_curve0,c_curve0,lambda_curve0,f_curve0_up,c_curve0_up,lambda_curve0_up,\
    f_curve0_low,c_curve0_low,lambda_curve0_low