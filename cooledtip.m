clear all 
close all


syms  w k r C1 C2 ua mueff P A B PI
 
% assume spherically symmetric
% u = R(r) / r   ==> laplacian( u ) =  R_rr(r) / r
% 
% steady state ==> w ( u - ua)  - k laplacian ( u )  = q
%
% where the laser source is given by
q = PI/4*P*mueff^2* exp (-mueff*r) / r

% solution is the homogenious plus a particular solution
% particular solution may be found from method of undetermined coefficients
% Rp = A *  exp (-mueff*r) + B  * r *  exp (-mueff*r) ;
A = PI/4*P*mueff^2/(w-k*mueff^2); 

% verify the particular solution
Rp = A *  exp (-mueff*r) ;
d2Rpdr2 = simple (diff(diff(Rp,r),r)); 
residual = simple(w * Rp   - k * d2Rpdr2 - q * r)


% general solution  to
%    w ( R(r)/r - ua)  - k R_rr(r)/r   = PI/4*P*mueff^2* exp (-mueff*r) / r
% of the form
up = Rp/r + ua;
uh = C1/r * exp( sqrt(w/k) * r ) + C2/r * exp( -sqrt(w/k) * r );
u = up + uh;
dudr =  diff(u,r);
% verify
rhs =  w*(u - ua) - k * diff(r^2 * dudr,r) / r^2 ;
simple(rhs-q)


% apply BC dirichlet at R1     no flux at R2
syms R1 R2 u0 
dupdr =  diff(up,r);
duhdr =  diff(uh,r);

%
%    | 1/r * exp( sqrt(w/k) * r )   1/r * exp( -sqrt(w/k) * r )   |
%A = |                                                            |
%    |     collect(duhdr,C1)               collect(duhdr,C2)      |
%

A = [ 1/R1 * exp( sqrt(w/k) * R1 ), 1/R1 * exp( -sqrt(w/k) * R1 ); (-1/R2^2*exp((w/k)^(1/2)*R2)+1/R2*(w/k)^(1/2)*exp((w/k)^(1/2)*R2)), (-1/R2^2*exp(-(w/k)^(1/2)*R2)-1/R2*(w/k)^(1/2)*exp(-(w/k)^(1/2)*R2))]

% dirichlet at R1     no flux at R2
%    | u0 -  up    |
%b = |             |
%    | 0 - dupdr   |
%
b = [ u0 - (PI/4*P*mueff^2/(w-k*mueff^2)*exp(-mueff*R1)/R1+ua ) ; -(-PI/4*P*mueff^2/(w-k*mueff^2)*mueff*exp(-mueff*R2)/R2-PI/4*P*mueff^2/(w-k*mueff^2)*exp(-mueff*R2)/R2^2)]

x  = A \ b;
x  = simple(x);

u = up  + [ 1/r * exp( sqrt(w/k) * r ) , 1/r * exp( -sqrt(w/k) * r )] *x;

%verify pde
dudr =  diff(u,r);
rhs =  w*(u - ua) - k * diff(r^2 * dudr,r) / r^2 ;
residual = simple(rhs-q);
%verify dirichlet BC
r  = 1 ;
R1 = 1 ;
checkdirichlet = simple(eval(u));
%verify neumann BC
r  = 5 ;
R2 = 5 ;
checkneumann = simple(eval(dudr));


ccode(u)

%>>s1 = 3.0/4.0/PI*P*mua*mutr/(w-k*mueff*mueff)*exp(-mueff*r)/r+ua;      s2 = s1;      s5 = 1/r*exp(sqrt(w/k)*r)*(-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*PI*R1*w+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*u0*PI*R1*k*mueff*mueff+3.0*sqrt(w/k)*R2*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)+4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*PI*R1*w-4.0*sqrt(w/k)*R2*exp(-sqrt(w/k)*R2)*ua*PI*R1*k*mueff*mueff-3.0*P*mua*mutr*mueff*R2*exp(-mueff*R2-sqrt(w/k)*R1)-3.0*P*mua*mutr*exp(-mueff*R2-sqrt(w/k)*R1)+4.0*exp(-sqrt(w/k)*R2)*ua*PI*R1*w-4.0*exp(-sqrt(w/k)*R2)*u0*PI*R1*w+4.0*exp(-sqrt(w/k)*R2)*u0*PI*R1*k*mueff*mueff+3.0*P*mua*mutr*exp(-sqrt(w/k)*R2-mueff*R1)-4.0*exp(-sqrt(w/k)*R2)*ua*PI*R1*k*mueff*mueff)/4.0;      s6 = exp(-sqrt(w/k)*(-R1+R2))/(-w+k*mueff*mueff)/PI/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      s4 = s5*s6;      s6 = 1/r*exp(-sqrt(w/k)*r)*exp(-sqrt(w/k)*(-R1+R2))/4.0;      s9 = 4.0*exp(sqrt(w/k)*R2)*u0*PI*R1*w-4.0*exp(sqrt(w/k)*R2)*u0*PI*R1*w*sqrt(w/k)*R2-4.0*exp(sqrt(w/k)*R2)*u0*PI*R1*k*mueff*mueff+4.0*exp(sqrt(w/k)*R2)*u0*PI*R1*k*mueff*mueff*sqrt(w/k)*R2-3.0*P*mua*mutr*exp(sqrt(w/k)*R2-mueff*R1)+3.0*P*mua*mutr*sqrt(w/k)*R2*exp(sqrt(w/k)*R2-mueff*R1)-4.0*exp(sqrt(w/k)*R2)*ua*PI*R1*w+4.0*exp(sqrt(w/k)*R2)*ua*PI*R1*w*sqrt(w/k)*R2+4.0*exp(sqrt(w/k)*R2)*ua*PI*R1*k*mueff*mueff-4.0*exp(sqrt(w/k)*R2)*ua*PI*R1*k*mueff*mueff*sqrt(w/k)*R2+3.0*P*mua*mutr*mueff*R2*exp(sqrt(w/k)*R1-mueff*R2)+3.0*P*mua*mutr*exp(sqrt(w/k)*R1-mueff*R2);      s10 = 1/(exp(-2.0*sqrt(w/k)*(-R1+R2))+sqrt(w/k)*R2*exp(-2.0*sqrt(w/k)*(-R1+R2))-1.0+sqrt(w/k)*R2);      s8 = s9*s10;      s9 = 1/PI/(-w+k*mueff*mueff);      s7 = s8*s9;      s5 = s6*s7;      s3 = s4+s5;      t0 = s2+s3;

% similarly for constant coeff
%syms q
%u = (q/w + ua) + C1/r * exp( sqrt(w/k) * r ) + C2/r * exp( -sqrt(w/k) * r );
%dudr =  diff(u,r);
%rhs =  w*(u - ua) - k * diff(r^2 * dudr,r) / r^2 ;
%simple(rhs-q)



