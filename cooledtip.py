from math import pi,sqrt,exp
import os

# all units should be meter, kelvin, kg, second
#
# cblood * omega * (u-ua) = [J / kg / K] * [kg /s /m^3 ] * K  = W / m^3
#
# k * d^2u/dx^2 = [W / m / K] * [ K /m^2 ] =  W / m^3
#
mua   = 0.45e+2
mus   = 47.0e+2
anfact= .9
mutr  = mua+mus*(1.0-anfact)
mueff = sqrt(3.0*mua*mutr)
ua    =  310
R1    =  .001
R2    =  .03
cblood = 3840.0

nstep = 100
def compradius(i):
   return R1 *(nstep - i ) / nstep + R2 * i / nstep 

def u(u0,k,w,P,r):
  s1 = 3.0/4.0/pi*P*mua*mutr/(w-k*mueff*mueff)*exp(-mueff*r)/r+ua;
  s2 = s1;
  s5 = 1/r*exp(sqrt(w/k)*r)*(-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*sqrt(w/k)*R2*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-3.0*P*mua*mutr*mueff*R2*exp(-mueff*R2-sqrt(w/k)*R1)-3.0*P*mua*mutr*exp(-mueff*R2-sqrt(w/k)*R1)+4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*w-4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*w+4.0*exp(-sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+3.0*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)-4.0*exp(-sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff)/4.0;
  s6 = exp(-sqrt(w/k)*(-R1+R2))/(-w+k*mueff*mueff)/pi/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);
  s4 = s5*s6;
  s6 = 1/r*exp(-sqrt(w/k)*r)*exp(-sqrt(w/k)*(-R1+R2))/4.0;
  s9 = 4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*w*sqrt(w/k)*R2-4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff+4.0*exp(sqrt(w/k)*R2)*u0*pi*R1*k*mueff*mueff*sqrt(w/k)*R2-3.0*P*mua*mutr*exp(sqrt(w/k)*R2-mueff*R1)+3.0*P*mua*mutr*sqrt(w/k)*R2*exp(sqrt(w/k)*R2-mueff*R1)-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*w*sqrt(w/k)*R2+4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff-4.0*exp(sqrt(w/k)*R2)*ua*pi*R1*k*mueff*mueff*sqrt(w/k)*R2+3.0*P*mua*mutr*mueff*R2*exp(sqrt(w/k)*R1-mueff*R2)+3.0*P*mua*mutr*exp(sqrt(w/k)*R1-mueff*R2);
  s10 = 1/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);
  s8 = s9*s10;
  s9 = 1/pi/(-w+k*mueff*mueff);
  s7 = s8*s9;
  s5 = s6*s7;
  s3 = s4+s5;
  t0 = s2+s3;
  return t0

u0range=[295,ua]# temp of applicator and baseline body temp

#Note that an acceptable range of thermal conductivity from CRC handbook
#is roughly .2-.6 [W/m/k]  . 5 is the average value for liver
krange = [.5 ]
#
#omega [kg/s/m^3] = .22 [ml/min/g] * 1 [g/cm^3]  
#                 = .22 [ml/min/cm^3]
#                 = .22 [g/min/cm^3]     (1g = 1ml for water)
#                 = .22 [g/min/cm^3] * [1kg/1000g] * [1min/60s] * [100cm/1m]^3
#                 = 3.6 [kg/s/m^3]
wrange = [3.6,36.]
#Prange = [0.0 ,5.0  ,15.0 ]
Prange = [5.0,20.]

paramlist =[ (u0,P_0, k_0,w_0,)
                    for u0      in u0range
                    for P_0     in  Prange
                    for k_0     in  krange
                    for w_0     in  wrange ]
RadRange = map( compradius, range(nstep+1))
datafile=open("cooledtip.dat"  ,"w")
for rad in RadRange:
  datafile.write("%f  "  % rad )
  for (u0,P_0, k_0,w_0) in paramlist:
    datafile.write("%f  "  % u(u0,k_0,cblood*w_0,P_0,rad)  )
  datafile.write("\n" )
datafile.close; datafile.flush()

gnuplotfile=open("cooledtip.plt"  ,"w")
gnuplotfile.write("""
# Gnuplot script file for plotting data in file "bloodperf.dat" 
#set term pslatex color size 7.5,5.25
#set   output "cooledtip.tex"
set   autoscale    # scale axes automatically   
unset log          # remove any log-scaling     
unset label        # remove any previous labels 
set xtic auto      # set xtics automatically    
set ytic auto      # set ytics automatically    
#set title "Nonlinear Blood cooledtip " 
set ylabel "temperature [degC]" 
set xlabel "radial distance [cm]" 
set key right  
q    = 1.0      
k1    = 0.01     
k2    = 0.10     
k3    = 1.00     
u43  = 43.0
u50  = 50.0
u0   = %f
ua   = %f
L    = 0.05       
#set xr [0:L]     
#set yr [-5.000000:105.000000]     
#f(x) = -q/k1*x*x/2 + ( (uL -u0)/L + q*L/2/k1)*x + u0
#g(x) = -q/k2*x*x/2 + ( (uL -u0)/L + q*L/2/k2)*x + u0
#h(x) = -q/k3*x*x/2 + ( (uL -u0)/L + q*L/2/k3)*x + u0
#set title "$ u(x) = -\\frac{q x^2}{2 k} + \\left( \\frac{u_L -u_0}{L} + \\frac{q L}{2k} \\right) x + u_0 $"
#set ytics  (  \
#"$u_0$"  u0   ,  \
#"$u_L$"  uL     )
#set xtics  (  \
#"0"  0   ,  \
#"$\\frac{L}{2} $"  L/2   ,  \
#"$L $"  L     )
plot \\
u50 notitle lt 0,\\
u43 notitle lt 0,\\
ua notitle lt 0,\\
u0 notitle lt 0,\\
"""  % (u0-273,ua-273) )
i = 1
nsize = len(paramlist)
for (u0,P_0, k_0,w_0) in paramlist:
    i=i+1
    gnuplotfile.write('"cooledtip.dat" using ($1*1000):($%d -273)  title "u0=%5.2f P_0= %5.2f k_0= %5.2f w_0= %5.2f"  w l lc %d lw 2 ' % (i,u0,P_0,k_0,w_0,i-1) )
    if (i != nsize +1):
       gnuplotfile.write(',\\\n' )
    else :
       gnuplotfile.write('\n' )

gnuplotfile.write("""
#unset output #to flush the output buffer
## create dvi file (dvi viewers DO NOT DISPLAY ROTATED TEXT)
#!latex -jobname cooledtip                                           \
#\\documentclass\{article\}                                           \
#\\usepackage\{amssymb,amsfonts,amsmath\}                             \
#\\usepackage\{nopageno\}                                             \
#\\usepackage[left=.5in,right=.5in,top=.5in,bottom=.5in]\{geometry\}  \
#\\begin\{document\}                                                  \
#\\LARGE                                                              \
#\\begin\{center\}                                                    \
#\\input\{cooledtip\}                                                \
#\\end\{center\}                                                      \
#\\end\{document\}
## convert to pdf (TEXT SHOULD BE ROTATED IN THE PDF WHERE APPROPRIATE)
#! dvipdf            cooledtip.dvi
## crop 
#! pdfcrop           cooledtip.pdf
#! mv                cooledtip-crop.pdf cooledtip.pdf
""")
gnuplotfile.close; gnuplotfile.flush()


os.system("gnuplot -persist cooledtip.plt -")
