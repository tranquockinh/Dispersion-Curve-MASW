from Paras import *
def MASW_plot_data(u, N, dx, L, T, Tmax, du, FigFontSize):
    SignalToPlot = np.zeros((np.size(u,0), N), dtype='float')
    fig, ax = plt.subplots(figsize=(5,6))
    
    for j in range(N):
        SignalToPlot[:,j] = u[:,j] + j * du * dx     
        
        plt.plot(SignalToPlot[:,j], T, 'k')
        
    ax.set_xlim(-5*dx*du,L*du+5*dx*du)
    ax.set_ylim(0,Tmax)
    ax.set_xlabel('Distance from source [m]', fontsize = FigFontSize)
    ax.set_ylabel('Time [s]', fontsize = FigFontSize)
    ax.set_xticks(np.arange(0,L*du+6*dx*du,6*dx*du))
    plt.gca().invert_yaxis()
    #plt.show()
    return SignalToPlot