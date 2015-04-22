from math import pi,sqrt,exp
import os

# all units should be meter, kelvin, kg, second
#
# cblood * omega * (u-ua) = [J / kg / K] * [kg /s /m^3 ] * K  = W / m^3
#
# k * d^2u/dx^2 = [W / m / K] * [ K /m^2 ] =  W / m^3
#

# use water viscosity and compressibility
Viscosity       = 8.9e-4  # \frac{kg}{m\; s}
Compressibility = 4.6e-10 #  \frac{m \; s^2}{kg}    
# guesstimate of permeability
PermeabilityList = [1.e-15,2.e-15,1.e-14] # m^2
PermeabilityList = [1.e-15,2.e-15] # m^2
DiffusivityList =  [ Perm/Viscosity/Compressibility  for Perm in PermeabilityList ] # m^2/s

Enthalpy = 45 # kJ/M

R1    =  .006
R2    =  .03
cblood = 3840.0
nstep = 100
def compradius(i):
   return R1 *(nstep - i ) / nstep + R2 * i / nstep 

def concentration(r,InjectionRate,InjectionLength,Diffusivity,c0):
  t0 = (-InjectionRate*(InjectionLength*InjectionLength)*R1*exp(-R2/InjectionLength)+InjectionRate*(InjectionLength*InjectionLength)*R1*exp(-r/InjectionLength)-InjectionRate*(InjectionLength*InjectionLength)*r*exp(-R1/InjectionLength)+InjectionRate*(InjectionLength*InjectionLength)*r*exp(-R2/InjectionLength)+Diffusivity*R1*c0*r-InjectionRate*InjectionLength*R1*R2*exp(-R2/InjectionLength)+InjectionRate*InjectionLength*R2*r*exp(-R2/InjectionLength))/(Diffusivity*R1*r);
  return t0

def temperature(r,InjectionRate,InjectionLength,Diffusivity,c0,u0,Conduction,Perfusion,ua):
  t0 = ua+((Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*exp(-r/InjectionLength))/(Diffusivity*(Perfusion-Conduction*1.0/(InjectionLength*InjectionLength)))-(Enthalpy*InjectionRate*InjectionLength*exp(-R2/InjectionLength)*(InjectionLength+R2))/(Diffusivity*Perfusion))/r+(Enthalpy*exp(-R1/InjectionLength)*exp(-R2/InjectionLength)*(InjectionRate*(InjectionLength*InjectionLength)*exp(R1/InjectionLength)-InjectionRate*(InjectionLength*InjectionLength)*exp(R2/InjectionLength)+InjectionRate*InjectionLength*R2*exp(R1/InjectionLength)+Diffusivity*R1*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)))/(Diffusivity*Perfusion*R1)-(1.0/(R1*R1)*exp(-R1/InjectionLength)*exp(-R2/InjectionLength)*exp(r*sqrt(Perfusion/Conduction))*(Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R1*R1)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R2*R2)*exp(R1/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2*R2)*exp(R1/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2)*exp(R1/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))+Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))+Conduction*Enthalpy*InjectionRate*InjectionLength*(R2*R2*R2)*exp(R1/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))+Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))-Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength)*Perfusion*R1*(R2*R2)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Conduction*Diffusivity*Enthalpy*(R1*R1*R1)*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))-Conduction*Diffusivity*Perfusion*(R1*R1*R1)*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Conduction*Diffusivity*Perfusion*(R1*R1*R1)*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R1*R1)*R2*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Diffusivity*Enthalpy*(InjectionLength*InjectionLength)*Perfusion*(R1*R1*R1)*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*R2*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)-Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*R2*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)-Conduction*Diffusivity*Enthalpy*(R1*R1*R1)*R2*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)-Conduction*Diffusivity*Perfusion*(R1*R1*R1)*R2*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Conduction*Diffusivity*Perfusion*(R1*R1*R1)*R2*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Diffusivity*Enthalpy*(InjectionLength*InjectionLength)*Perfusion*(R1*R1*R1)*R2*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)))/(Diffusivity*Perfusion*r*(Conduction-(InjectionLength*InjectionLength)*Perfusion)*(exp(R1*sqrt(Perfusion/Conduction)*2.0)-exp(R2*sqrt(Perfusion/Conduction)*2.0)+R2*exp(R1*sqrt(Perfusion/Conduction)*2.0)*sqrt(Perfusion/Conduction)+R2*exp(R2*sqrt(Perfusion/Conduction)*2.0)*sqrt(Perfusion/Conduction)))+(1.0/(R1*R1)*exp(-R1/InjectionLength)*exp(-R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))*exp(R2*sqrt(Perfusion/Conduction))*exp(-r*sqrt(Perfusion/Conduction))*(Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R2*R2)*exp(R1/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R1*R1)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2*R2)*exp(R1/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))-Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2)*exp(R1/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength*InjectionLength)*Perfusion*(R2*R2)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Conduction*Enthalpy*InjectionRate*InjectionLength*(R2*R2*R2)*exp(R1/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))+Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))+Enthalpy*InjectionRate*(InjectionLength*InjectionLength*InjectionLength)*Perfusion*R1*(R2*R2)*exp(R2/InjectionLength)*exp(R1*sqrt(Perfusion/Conduction))-Conduction*Diffusivity*Enthalpy*(R1*R1*R1)*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Conduction*Diffusivity*Perfusion*(R1*R1*R1)*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))+Conduction*Diffusivity*Perfusion*(R1*R1*R1)*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Conduction*Enthalpy*InjectionRate*(InjectionLength*InjectionLength)*(R1*R1)*R2*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Diffusivity*Enthalpy*(InjectionLength*InjectionLength)*Perfusion*(R1*R1*R1)*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))-Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*R2*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Diffusivity*(InjectionLength*InjectionLength)*(Perfusion*Perfusion)*(R1*R1*R1)*R2*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Conduction*Diffusivity*Enthalpy*(R1*R1*R1)*R2*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)+Conduction*Diffusivity*Perfusion*(R1*R1*R1)*R2*u0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)-Conduction*Diffusivity*Perfusion*(R1*R1*R1)*R2*ua*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)-Diffusivity*Enthalpy*(InjectionLength*InjectionLength)*Perfusion*(R1*R1*R1)*R2*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)*exp(R2*sqrt(Perfusion/Conduction))*sqrt(Perfusion/Conduction)))/(Diffusivity*Perfusion*r*(Conduction-(InjectionLength*InjectionLength)*Perfusion)*(exp(R1*sqrt(Perfusion/Conduction)*2.0)-exp(R2*sqrt(Perfusion/Conduction)*2.0)+R2*exp(R1*sqrt(Perfusion/Conduction)*2.0)*sqrt(Perfusion/Conduction)+R2*exp(R2*sqrt(Perfusion/Conduction)*2.0)*sqrt(Perfusion/Conduction)));
  return t0

u0range=[100.5]# temp of applicator and baseline body temp

#Note that an acceptable range of thermal conductivity from CRC handbook
#is roughly .2-.6 [W/m/k]  . 5 is the average value for liver
krange = [.5 ]
#
#omega [kg/s/m^3] = .22 [ml/min/g] * 1 [g/cm^3]  
#                 = .22 [ml/min/cm^3]
#                 = .22 [g/min/cm^3]     (1g = 1ml for water)
#                 = .22 [g/min/cm^3] * [1kg/1000g] * [1min/60s] * [100cm/1m]^3
#                 = 3.6 [kg/s/m^3]
cblood = 3840
wrange = [ 6.0]
PerfusionList = [ w_0*cblood for w_0 in wrange]
arterytemp=37.0 # degC

massinjectionList =  [  200.]  # kg/m^2/s
lengthinjectionList = [ .010]  # m
BoundaryConcentration = 7.5 # Molar = mol/L
MolarMass = 82.0343 # (kg L)/(m^3 mol)
y2tics = [ '"%4.1f" %f' % (molar,molar*MolarMass ) for  molar in [1,3,5,7.5,9]]
boundarydensity     = BoundaryConcentration *MolarMass # kg/m^3
paramlist =[ (u0,massinjection, lengthinjection ,diffcoef , k_0,perf)
                  for u0              in u0range
                  for massinjection   in massinjectionList 
                  for lengthinjection in lengthinjectionList 
                  for diffcoef        in DiffusivityList 
                  for k_0             in krange
                  for perf            in PerfusionList ]
RadRange = map( compradius, range(nstep+1))
datafile=open("ficksdiffusion.dat"  ,"w")
for rad in RadRange:
  datafile.write("%f  "  % rad )
  for (u0 ,massinjection,lengthinjection,diffcoef , k_0,perf) in paramlist:
    cvalue = concentration(rad,massinjection,lengthinjection,diffcoef,boundarydensity ) 
    uvalue = temperature(rad,massinjection,lengthinjection,diffcoef,boundarydensity,u0,k_0,perf,arterytemp) 
    datafile.write("%f  %f "  %  ( cvalue , uvalue ) ) 
  datafile.write("\n" )
datafile.close; datafile.flush()

gnuplotfile=open("ficksdiffusion.plt"  ,"w")
gnuplotfile.write("""
# Gnuplot script file for plotting data in file "bloodperf.dat" 
set term pslatex color size 7.5,5.25
set   output "ficksdiffusion.tex"
set   autoscale    # scale axes automatically   
unset log          # remove any log-scaling     
unset label        # remove any previous labels 
set xtic auto      # set xtics automatically    
set ytic auto  nomirror     # set ytics automatically    
set y2tic (%s)  nomirror    # set y2tics automatically    
set ylabel "Temperature [$^{\\\\circ}$C]" 
set y2label "Concentration [M]" 
set xlabel "radial distance [mm]" 
set key right  spacing 1.5
q    = 1.0      
k1    = 0.01     
k2    = 0.10     
k3    = 1.00     
u43  = 43.0
u50  = 50.0
u0   = %f
ua   = %f
L    = 0.05       
set xr [0:30]     
set label "$ D \\\\Delta c = c0 \\\\frac{\\\\exp{d\\\\; r}}{r}$ " at 20,60 
#set yr [-5.000000:105.000000]     
set y2r [0.000000:700.000000]     
#f(x) = -q/k1*x*x/2 + ( (uL -u0)/L + q*L/2/k1)*x + u0
#g(x) = -q/k2*x*x/2 + ( (uL -u0)/L + q*L/2/k2)*x + u0
#h(x) = -q/k3*x*x/2 + ( (uL -u0)/L + q*L/2/k3)*x + u0
#set title "$ u(x) = -\\\\frac{q x^2}{2 k} + \\\\left( \\\\frac{u_L -u_0}{L} + \\\\frac{q L}{2k} \\\\right) x + u_0 $"
plot \\
u50 notitle lt 0,\\
u43 notitle lt 0,\\
ua notitle lt 0,\\
u0 notitle lt 0,\\
"""  % (",".join(y2tics),u0,arterytemp) )
i = 1
nsize = len(paramlist)
for (u0,massinjection,lengthinjection,diffcoef , k_0,perf) in paramlist:
    i=i+1
    gnuplotfile.write('"ficksdiffusion.dat" using ($1*1000):($%d )  title "$D$= %10.3e$m^2/s$ $k$= %10.3e "  w l lc %d lw 2 lt 1  axes x1y2' % ( 2*(i-1)  ,diffcoef,k_0,i-1) )
    gnuplotfile.write(',\\\n' )
    gnuplotfile.write('"ficksdiffusion.dat" using ($1*1000):($%d )  notitle                                  w l lc %d lw 2 lt 2  axes x1y1' % ( 2*(i-1)+1,i-1) )
    if (i != nsize +1):
       gnuplotfile.write(',\\\n' )
    else :
       gnuplotfile.write('\n' )

gnuplotfile.write("""
unset output #to flush the output buffer
# create dvi file (dvi viewers DO NOT DISPLAY ROTATED TEXT)
!latex -jobname ficksdiffusion                                       \
\\\\documentclass{article}                                             \
\\\\usepackage{amssymb,amsfonts,amsmath}                               \
\\\\usepackage{nopageno}                                               \
\\\\usepackage[left=.5in,right=.5in,top=1.0in,bottom=.5in]{geometry}   \
\\\\begin{document}                                                    \
\\\\LARGE                                                              \
\\\\begin{center}                                                      \
\\\\input{ficksdiffusion}                                              \
\\\\end{center}                                                        \
\\\\end{document}                                                      \
# convert to pdf (TEXT SHOULD BE ROTATED IN THE PDF WHERE APPROPRIATE)
! dvipdf            ficksdiffusion.dvi
# crop 
! pdfcrop           ficksdiffusion.pdf
! mv                ficksdiffusion-crop.pdf ficksdiffusion.pdf
""")
gnuplotfile.close; gnuplotfile.flush()


os.system("gnuplot ficksdiffusion.plt ")
#os.system("gnuplot -persist ficksdiffusion.plt -")
