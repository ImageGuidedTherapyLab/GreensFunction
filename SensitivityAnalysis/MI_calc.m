close all;
clear;

%% Parameters
ndim=240;
%k=[0,0];
mu_hat=200;
sigma_nu=1;
sigma_mu=0.001;

M(1:ndim,1:ndim)=1;
%T2star=[10,70,100,80];
T2star(1:ndim,1:ndim) = 0.1;
delta_omega0=0;
T_E=0.025;
gamma=42.58;
alpha=-0.0102;
B_0=1.5;

Power=12;

r_1=0.0075;
r_2=1;

c_p=3640;
c_blood=3840;
rho=1045;
k=0.5270;
omega=6;%9;

%% Calculate Boundary Condition Constants
[C_1,C_2,const_params,const_params2] = sens_const(omega,c_blood,k,mu_hat,Power,r_1,r_2);

%%

rVar=zeros(240);
for i=1:240
    for j=1:240
        rVar(i,j)=sqrt((i-120.5)^2+(j-120.5)^2)/10000;
    end
end

%% Calculate temperature field
create_temp_field_DM
load('all_opt_mu_small.mat');
u1 = all_opt_fig(:,:,1);
u2 = all_opt_fig(:,:,2);
%u=temp_DM(omega,c_blood,k,mu_hat,Power,r_1,r_2,rVar);

%% Calculate temperature sensitivity
%du_dmu=temp_sens(omega,c_blood,k,mu_hat,Power,r_1,r_2,rVar);
%du_dmu=sens_soln(C_1,C_2,const_params,const_params2,omega,c_blood,k,mu_hat,r_1,rVar);
du1_dmu=all_opt_fig_s(:,:,1);
du2_dmu=all_opt_fig_s(:,:,2);

%% Evaluate signal model, G(mu,k)
s = T_E./T2star + 1i*(2*pi*gamma*alpha*B_0*T_E*u1 + T_E*delta_omega0);
ds_dmu = 1i*2*pi*gamma*alpha*B_0*T_E*du1_dmu;

z_pretrans = M.*exp(-s);
z = fftshift(fft2(z_pretrans));
b = abs(z);

dz_dmu_pretrans = z_pretrans.*ds_dmu;
z_prime = fftshift(fft2(dz_dmu_pretrans));
%a = z.*z_prime./b;
%a = real(z.*z_prime)./b;
a = (real(z).*real(z_prime) + imag(z).*imag(z_prime))./b;

%% Calculate mutual information from analytical solution involving normal distributions
MI = log(sqrt(sigma_nu^2 + a.^2*sigma_mu^2)/sigma_nu);

%% Calculate approximate sensitivity for verification
delu_delmu = u2 - u1;

s = T_E./T2star + 1i*(2*pi*gamma*alpha*B_0*T_E*u2 + T_E*delta_omega0);
ds_dmu = 1i*2*pi*gamma*alpha*B_0*T_E*delu_delmu;

z_pretrans = M.*exp(-s);
z = fftshift(fft2(z_pretrans));
b = abs(z);

dz_dmu_pretrans = z_pretrans.*ds_dmu;
z_prime = fftshift(fft2(dz_dmu_pretrans));
%a = z.*z_prime./b;
%a = real(z.*z_prime)./b;
a = (real(z).*real(z_prime) + imag(z).*imag(z_prime))./b;

MI_approx = log(sqrt(sigma_nu^2 + a.^2*sigma_mu^2)/sigma_nu);

%% Display results
figure; imagesc(u1);
figure; imagesc(du1_dmu);
figure; imagesc(delu_delmu);
figure; imagesc(MI);
figure; imagesc(MI_approx);