% This function calculates du/dmu, the temperature sensitivity to the optical
% attenuation coefficient, using the analytical solution to the steady
% state of the Pennes Bioheat Equation, which is a non-homogeneous 
% second-order ODE. 

function [sensitivity] = temp_sens(omega,c_blood,k,mu_eff,P,r_1,r_2,r)

%% Calculate constants in homogeneous solution
const_params = sqrt(omega*c_blood/k);
A_plus = exp(const_params*r_1)/r_1;
A_minus = exp(-const_params*r_1)/r_1;
B_plus = (const_params*r_2*exp(const_params*r_2) - exp(const_params*r_2))/r_2^2;
B_minus = (-const_params*r_2*exp(-const_params*r_2) - exp(-const_params*r_2))/r_2^2;

const_params2 = mu_eff*P/(4.0*pi*(omega*c_blood - k*mu_eff^2)^2);
D_1 = const_params2*(k*mu_eff^3*r_1 + c_blood*(2.0 - mu_eff*r_1)*omega)*exp(-mu_eff*r_1);
D_2 = const_params2*((k*mu_eff^3 - c_blood*mu_eff*omega)*exp(-mu_eff*r_2) ...
        - mu_eff*(k*mu_eff^3*r_2 + c_blood*(2.0 - mu_eff*r_2)*omega)*exp(-mu_eff*r_2));

C_2 = (A_plus*D_2/B_plus - D_1)/(B_minus*A_plus/B_plus - A_minus);
C_1 = (D_1 - A_minus*C_2)/A_plus;

%% Homogeneous solution
soln_h = C_1*exp(const_params*r)./r + C_2*exp(-const_params*r)./r;

%% Particular solution
soln_p = const_params2*(k*mu_eff^3*r + c_blood*(2.0 - mu_eff*r)*omega)*exp(-mu_eff*r)./r;

%% General solution
sensitivity = soln_h + soln_p;

%% Boundary conditions
sensitivity.*(r>=r_1);