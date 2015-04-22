clear all 
close all

syms Diffusivity InjectionRate InjectionLength A r C1 C2 ; 

%Diffusivity = PhaseMobility/Compressibility = K/mu/Compressibility

% assume spherically symmetric
% c = R(r) / r   ==> laplacian( c ) =  R_rr(r) / r
% 
% steady state ==>  - Diffusivity laplacian ( c ) + MassSource = 0
%
disp(' mass source  [kg/m^3/s] ');
MassSource = InjectionRate / r  * exp (-r/InjectionLength) 

% solution is the homogenious plus a particular solution
% particular solution may be found from method of undetermined coefficients
% cRp = A *  exp (-mueff*r) + B  * r *  exp (-mueff*r) ;
A = ( InjectionRate* InjectionLength^2 )/Diffusivity  

disp(' verify the particular solution ');
cRp = A *  exp (-r/InjectionLength) ;
d2Rpdr2 = simplify(diff(diff(cRp,r),r)); 
residual = simplify( - Diffusivity* d2Rpdr2 + MassSource * r);




% general solution  to
%     - K / mu / Compressibility  R_rr(r)/r  +  InjectionRate/r   * exp (-r/InjectionLength)  = 0 
% of the form
cp = cRp/r 
ch = C1 + C2/r;
c = cp + ch;
dcdr =  diff(c,r);
% verify
rhs =   - Diffusivity * diff(r^2 * dcdr,r) / r^2 ;
concentrationresidual = simplify(rhs+MassSource );


% apply BC dirichlet at R1     no flux at R2
syms R1 R2 c0 
dcpdr =  diff(cp,r);
dchdr =  diff(ch,r);

%
%          |           1          1/R1              |
% matrix = |                                        |
%          |           0         -1/R2^2            |
%

matrix = [ 1, 1/R1 ; 0 , -1/R2^2 ] 

% dirichlet at R1     no flux at R2
%    | c0 -  cp    |
%b = |             |
%    | 0 - dcpdr   |
%
b = [ c0 - ( InjectionRate* InjectionLength^2 )/Diffusivity *  exp(-R1/InjectionLength)/R1 ; (InjectionRate*InjectionLength*exp(-R2/InjectionLength))/(Diffusivity*R2) + (InjectionRate*InjectionLength^2*exp(-R2/InjectionLength))/(Diffusivity*R2^2) ]

x  = matrix \ b;

c = simplify(cp  + [ 1 , 1/r ] *x);

%verify pde
dcdr =  diff(c,r);
rhs =  Diffusivity * diff(r^2 * dcdr,r) / r^2 ;
residual = simplify(rhs+MassSource);
%% %verify dirichlet BC
%% r  = 1 ;
%% R1 = 1 ;
%% checkdirichlet = simplify(eval(c));
%% %verify neumann BC
%% r  = 5 ;
%% R2 = 5 ;
%% checkneumann = simplify(eval(dcdr));

%% 
ccode(c)
%>> t0 = (-InjectionRate*(InjectionLength*InjectionLength)*R1*exp(-R2/InjectionLength)+InjectionRate*(InjectionLength*InjectionLength)*R1*exp(-r/InjectionLength)-InjectionRate*(InjectionLength*InjectionLength)*r*exp(-R1/InjectionLength)+InjectionRate*(InjectionLength*InjectionLength)*r*exp(-R2/InjectionLength)+Diffusivity*R1*c0*r-InjectionRate*InjectionLength*R1*R2*exp(-R2/InjectionLength)+InjectionRate*InjectionLength*R2*r*exp(-R2/InjectionLength))/(Diffusivity*R1*r);


syms  Perfusion Conduction ua Enthalpy P B PI r R1 R2 C3 C4 u0
 
% assume spherically symmetric
% u = R(r) / r   ==> laplacian( u ) =  R_rr(r) / r
% 
% steady state ==> Perfusion ( u - ua)  - Conduction laplacian ( u )  = Enthalpy  * c
%
% where the laser source is given by
q = Enthalpy  * (x(2) /r +cp)
%q = Enthalpy  * ( cp)

% solution is the homogenious plus a particular solution
% particular solution may be found from method of undetermined coefficients
% uRp = A *  exp (-mueff*r) + B  * r *  exp (-mueff*r) ;

% verify the particular solution
B =  Enthalpy / (Perfusion-Conduction/InjectionLength^2) 
uRp =  Enthalpy* x(2)/Perfusion  + B *  cp *r   ;
%uRp =   B *  cp* r ;
d2Rpdr2  = simplify( diff(diff(uRp,r),r)); 
residual = simplify(Perfusion * uRp   - Conduction * d2Rpdr2 - q * r);


% general solution  to
%    Perfusion ( R(r)/r - ua)  - Conduction R_rr(r)/r   = 3/4/PI*P*mua*mutr * exp (-mueff*r) / r
% of the form
up = uRp/r + ua + Enthalpy/Perfusion *x(1);
uh = C3/r * exp( sqrt(Perfusion/Conduction) * r ) + C4/r * exp( -sqrt(Perfusion/Conduction) * r );
u = up + uh;
dudr  =  diff(u,r);
dupdr =  diff(up,r);
% verify
rhs =  Perfusion*(u - ua) - Conduction * diff(r^2 * dudr,r) / r^2 ;
residualu = simplify(rhs-q - Enthalpy*x(1));

%
%              | 1/r * exp( sqrt(Perfusion/Conduction) * r )   1/r * exp( -sqrt(Perfusion/Conduction) * r )   |
% tempmatrix = |                                                                                              |
%              |                 collect(duhdr,C1)                                  collect(duhdr,C2)         |
%

tempmatrix = [ 1/R1 * exp( sqrt(Perfusion/Conduction) * R1 ), 1/R1 * exp( -sqrt(Perfusion/Conduction) * R1 ); (-1/R2^2*exp((Perfusion/Conduction)^(1/2)*R2)+1/R2*(Perfusion/Conduction)^(1/2)*exp((Perfusion/Conduction)^(1/2)*R2)), (-1/R2^2*exp(-(Perfusion/Conduction)^(1/2)*R2)-1/R2*(Perfusion/Conduction)^(1/2)*exp(-(Perfusion/Conduction)^(1/2)*R2))]

% dirichlet at R1     no flux at R2
%    | u0 -  up    |
%b = |             |
%    | 0 - dupdr   |
%
temprhs = [ u0 - ua + ((Enthalpy*InjectionRate*InjectionLength^2*exp(-R1/InjectionLength))/(Diffusivity*(Perfusion - Conduction/InjectionLength^2)) - (Enthalpy*InjectionRate*InjectionLength*exp(-R2/InjectionLength)*(InjectionLength + R2))/(Diffusivity*Perfusion))/R1 + (Enthalpy*exp(-R1/InjectionLength)*exp(-R2/InjectionLength)*(InjectionRate*InjectionLength^2*exp(R1/InjectionLength) - InjectionRate*InjectionLength^2*exp(R2/InjectionLength) + InjectionRate*InjectionLength*R2*exp(R1/InjectionLength) + Diffusivity*R1*c0*exp(R1/InjectionLength)*exp(R2/InjectionLength)))/(Diffusivity*Perfusion*R1)  ; 
-( - ((Enthalpy*InjectionRate*InjectionLength^2*exp(-R1/InjectionLength))/(Diffusivity*(Perfusion - Conduction/InjectionLength^2)) - (Enthalpy*InjectionRate*InjectionLength*exp(-R2/InjectionLength)*(InjectionLength + R2))/(Diffusivity*Perfusion))/R1^2 - (Enthalpy*InjectionRate*InjectionLength*exp(-R1/InjectionLength))/(Diffusivity*R1*(Perfusion - Conduction/InjectionLength^2)) ) ]


y  = tempmatrix \ temprhs;


u = up  + [ 1/r * exp( sqrt(Perfusion/Conduction) * r ) , 1/r * exp( -sqrt(Perfusion/Conduction) * r )] *y;

ccode(u)
