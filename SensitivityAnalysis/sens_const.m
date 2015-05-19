%% Calculate the constants from boundary conditions for sensitivity solution

function [C_1,C_2,const_params,const_params2] = sens_const(omega,c_blood,k,mu_eff,P,r_1,r_2)

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