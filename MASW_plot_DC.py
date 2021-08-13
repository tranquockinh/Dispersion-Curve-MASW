from Paras import *
def MASW_plot_DC(f_curve0,c_curve0,lambda_curve0,f_curve0_up,c_curve0_up,\
    lambda_curve0_up,f_curve0_low,c_curve0_low,lambda_curve0_low):
    
    ## types of plots to display
    plot_f_c ='f_c'
    plot_lambda_c = 'lambda_c'
    
    
    if plot_f_c == 'f_c':
        fig, ax = plt.subplots()
        ax.plot(f_curve0,c_curve0,'ko-')
        ax.plot(f_curve0_up,c_curve0_up,'r+--')
        ax.plot(f_curve0_low,c_curve0_low,'r+--')
        ax.set_xlabel('Frequency (Hz)', fontsize = FigFontSize)
        ax.set_ylabel('Reyleigh wave velocity [m/s]', fontsize = FigFontSize)
        ax.legend(['Fundamental DC','Upper BC','Lower BC'],loc='upper right')
    if plot_lambda_c == 'lambda_c':
        
        fig, ax = plt.subplots()
        ax.plot(c_curve0,lambda_curve0, 'ko-')
        ax.plot(c_curve0_up,lambda_curve0_up,'r+--')
        ax.plot(c_curve0_low,lambda_curve0_low,'r+--')
        plt.gca().invert_yaxis()
        ax.set_xlabel('Reyleigh wave velocity [m/s]', fontsize = FigFontSize)
        ax.set_ylabel('Wavelength [m]', fontsize = FigFontSize)
        ax.legend(['Fundamental DC','Upper BC','Lower BC'],loc='upper right')
    plt.show()
