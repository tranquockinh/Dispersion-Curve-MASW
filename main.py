
from Paras import *
from MASW_plot_data import MASW_plot_data
from MASW_Despersion_image import MASW_Despersion_image
from MASW_plot_Dispersion_image import MASW_plot_Dispersion_image
from MASW_Extract_Dispersion_Curve import MASW_Extract_Dispersion_Curve
from MASW_plot_DC import MASW_plot_DC

SignalToPlot = MASW_plot_data(u,N,dx,L,T,Tmax,du,FigFontSize)

f,c,A = MASW_Despersion_image(u,N,fs,cT_min,cT_max,delta_cT)
#print('f\n{}\nc\n{}\nA\n{}'.format(f,c,A))
#print('f_shape\n{}\nc_shape\n{}\nA_shape\n{}'.format(f.shape,c.shape,A.shape))

# Plot dispersion image
valMin, valMax, IdxMin, IdxMax, fplot, cplot, Aplot  = MASW_plot_Dispersion_image(f,c,A)
print(fplot.shape, cplot.shape, Aplot.shape)
print('fplot:\n{}\ncplot:\n{}\nAplot:\n{}'.format(fplot, cplot, Aplot))

# Inversion analysis
            
f_curve0,c_curve0,lambda_curve0,f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,\
c_curve0_low,lambda_curve0_low = MASW_Extract_Dispersion_Curve(f,c,A) 

    
MASW_plot_DC(f_curve0,c_curve0,lambda_curve0,f_curve0_up,c_curve0_up,lambda_curve0_up,\
f_curve0_low,c_curve0_low,lambda_curve0_low)
