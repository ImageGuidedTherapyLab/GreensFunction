% This function calculates du/dmu, the temperature sensitivity to the optical
% attenuation coefficient, using the analytical solution to the steady
% state of the Pennes Bioheat Equation, which is a non-homogeneous 
% second-order ODE. 

function [sensitivity] = sens_soln(C_1,C_2,const_params,const_params2,omega,c_blood,k,mu_eff,r_1,r)

%% Homogeneous solution
soln_h = C_1*exp(const_params.*r)./r + C_2*exp(-const_params.*r)./r;

%% Particular solution
soln_p = const_params2*(k*mu_eff^3.*r + c_blood*(2.0 - mu_eff.*r)*omega)*exp(-mu_eff.*r)./r;

%% General solution
sensitivity = soln_h + soln_p;

%% Boundary conditions
sensitivity.*(r>=r_1);