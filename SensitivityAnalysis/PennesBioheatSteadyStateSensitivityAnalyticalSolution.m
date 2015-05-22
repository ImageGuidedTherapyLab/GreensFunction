%% This script generates a file in C which contains the analytic solution to the steady state sensitivity analysis of the Pennes Bioheat Equation.

syms sensitivity(rVar) kCond wPerf cblood mueff r1 r2;

Dsensitivity = diff(sensitivity);

sensitivity(rVar) = dsolve(kCond*diff(sensitivity,rVar,rVar) + 2.0*kCond*diff(sensitivity,rVar)/rVar - wPerf*cblood*sensitivity == mueff^2*exp(-mueff*rVar)/(4.0*pi) - mueff*exp(-mueff*rVar)/(2.0*pi*rVar), sensitivity(r1) == 0, Dsensitivity(r2) == 0);

sensitivity(rVar) = simplify(sensitivity);

ccode(sensitivity,'file','PennesBioheatSteadyStateSensitivityAnalyticalSolution');