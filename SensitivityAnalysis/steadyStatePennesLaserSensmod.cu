/*
 * Example Matlab cuda kernel interface.
 */


__device__
void pointSource(double rVar, double r1, double r2, double wPerf, double cblood, double kCond, double mueff, double u0, double ua, double Power, double *temperature )
{

   double pi = 3.141592653589793;
//*temperature  =ua-(exp(-mueff*rVar)*powf(mueff,2)*Power)/(4*kCond*powf(mueff,2)*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar);

*temperature=ua-(exp(-mueff*rVar)*powf(mueff,2)*Power)/(4*kCond*powf(mueff,2)*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar)-exp(mueff*rVar+r2*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(rVar,2)*(-1+r2*sqrt((cblood*wPerf)/kCond))+4*exp(mueff*(r1+rVar)+r2*sqrt((cblood*wPerf)/kCond))*pi*r1*powf(rVar,2)*(u0-ua)*(-kCond*powf(mueff,2)+cblood*wPerf)*(-1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*powf(rVar,3)*(-kCond*powf(mueff,2)+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))))+(exp(-mueff*(r1+rVar)+(2*r1+rVar)*sqrt((cblood*wPerf)/kCond))*(-exp(mueff*r1+r2*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar)-exp(mueff*rVar+r1*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(rVar,2)*(1+r2*sqrt((cblood*wPerf)/kCond))-4*exp(mueff*(r1+rVar)+r1*sqrt((cblood*wPerf)/kCond))*pi*r1*powf(rVar,2)*(u0-ua)*(kCond*powf(mueff,2)-cblood*wPerf)*(1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*powf(r1,2)*powf(rVar,3)*(-kCond*powf(mueff,2)+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))));

//	*temperature = ua-(exp(-mueff*rVar)*mueff*mueff*Power)/(4*kCond*mueff*mueff*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*mueff*mueff*Power*r2*r2*(1+mueff*rVar)-exp(mueff*rVar+r2*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*rVar*rVar*(-1+r2*sqrt((cblood*wPerf)/kCond))+4*exp(mueff*(r1+rVar)+r2*sqrt((cblood*wPerf)/kCond))*pi*r1*rVar*rVar*(u0-ua)*(-kCond*mueff*mueff+cblood*wPerf)*(-1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*rVar*rVar*rVar*(-kCond*mueff*mueff+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))))+(exp(-mueff*(r1+rVar)+(2*r1+rVar)*sqrt((cblood*wPerf)/kCond))*(-exp(mueff*r1+r2*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*r2*r2*(1+mueff*rVar)-exp(mueff*rVar+r1*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*rVar*rVar*(1+r2*sqrt((cblood*wPerf)/kCond))-4*exp(mueff*(r1+rVar)+r1*sqrt((cblood*wPerf)/kCond))*pi*r1*rVar*rVar*(u0-ua)*(kCond*mueff*mueff-cblood*wPerf)*(1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*r1*r1*rVar*rVar*rVar*(-kCond*mueff*mueff+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))));

//   *temperature = ua+(P*PI_Var*(mueff*mueff)*exp(-mueff*r)*(1.0/4.0))/(r*(w-k*(mueff*mueff)))-(exp(-R1*mueff-R2*mueff)*exp(r*sqrt(w/k))*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*R2*(mueff*mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-R1*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)))-(exp(R1*sqrt(w/k))*exp(R2*sqrt(w/k))*exp(-r*sqrt(w/k))*exp(-R1*mueff)*exp(-R2*mueff)*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)+P*PI_Var*R2*(mueff*mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)+R1*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)));
}
__device__
void pointSourceSensLong(double rVar, double r1, double r2, double wPerf, double cblood, double kCond, double mueff, double u0, double ua, double Power, double *sensitivity )
{
    double pi = 3.141592653589793;
// General solution
*sensitivity = ((((mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*(kCond*mueff*mueff*mueff*r1 + cblood*(2.0 - mueff*r1)*wPerf)*exp(-mueff*r1)) - (exp(-sqrt(wPerf*cblood/kCond)*r1)/r1)*(((exp(sqrt(wPerf*cblood/kCond)*r1)/r1)*((mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*((kCond*mueff*mueff*mueff - cblood*mueff*wPerf)*exp(-mueff*r2) - mueff*(kCond*mueff*mueff*mueff*r2 + cblood*(2.0 - mueff*r2)*wPerf)*exp(-mueff*r2)))/((sqrt(wPerf*cblood/kCond)*r2*exp(sqrt(wPerf*cblood/kCond)*r2) - exp(sqrt(wPerf*cblood/kCond)*r2))/(r2*r2)) - ((mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*(kCond*mueff*mueff*mueff*r1 + cblood*(2.0 - mueff*r1)*wPerf)*exp(-mueff*r1)))/(((-sqrt(wPerf*cblood/kCond)*r2*exp(-sqrt(wPerf*cblood/kCond)*r2) - exp(-sqrt(wPerf*cblood/kCond)*r2))/(r2*r2))*(exp(sqrt(wPerf*cblood/kCond)*r1)/r1)/((sqrt(wPerf*cblood/kCond)*r2*exp(sqrt(wPerf*cblood/kCond)*r2) - exp(sqrt(wPerf*cblood/kCond)*r2))/(r2*r2)) - (exp(-sqrt(wPerf*cblood/kCond)*r1)/r1))))/(exp(sqrt(wPerf*cblood/kCond)*r1)/r1))*exp(sqrt(wPerf*cblood/kCond)*rVar)/rVar + (((exp(sqrt(wPerf*cblood/kCond)*r1)/r1)*((mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*((kCond*mueff*mueff*mueff - cblood*mueff*wPerf)*exp(-mueff*r2) - mueff*(kCond*mueff*mueff*mueff*r2 + cblood*(2.0 - mueff*r2)*wPerf)*exp(-mueff*r2)))/((sqrt(wPerf*cblood/kCond)*r2*exp(sqrt(wPerf*cblood/kCond)*r2) - exp(sqrt(wPerf*cblood/kCond)*r2))/(r2*r2)) - ((mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*(kCond*mueff*mueff*mueff*r1 + cblood*(2.0 - mueff*r1)*wPerf)*exp(-mueff*r1)))/(((-sqrt(wPerf*cblood/kCond)*r2*exp(-sqrt(wPerf*cblood/kCond)*r2) - exp(-sqrt(wPerf*cblood/kCond)*r2))/(r2*r2))*(exp(sqrt(wPerf*cblood/kCond)*r1)/r1)/((sqrt(wPerf*cblood/kCond)*r2*exp(sqrt(wPerf*cblood/kCond)*r2) - exp(sqrt(wPerf*cblood/kCond)*r2))/(r2*r2)) - (exp(-sqrt(wPerf*cblood/kCond)*r1)/r1)))*exp(-sqrt(wPerf*cblood/kCond)*rVar)/rVar + (mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff)))*(kCond*mueff*mueff*mueff*rVar + cblood*(2.0 - mueff*rVar)*wPerf)*exp(-mueff*rVar)/rVar;
//printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n",mueff,Power,wPerf,cblood,kCond,r1,rVar,*sensitivity);
}
__device__
void pointSourceSens(double rVar, double r1, double r2, double wPerf, double cblood, double kCond, double mueff, double u0, double ua, double Power, double *sensitivity )
{

   double pi = 3.141592653589793;

// Calculate constants in homogeneous solution
double const_params;
double const_params2;
double A_plus;
double A_minus;
double B_plus;
double B_minus;
double C_1;
double C_2;
double D_1;
double D_2;
double soln_h;
double soln_p;
double sens;

const_params = sqrt(wPerf*cblood/kCond);
A_plus = exp(const_params*r1)/r1;
A_minus = exp(-const_params*r1)/r1;
B_plus = (const_params*r2*exp(const_params*r2) - exp(const_params*r2))/(r2*r2);
B_minus = (-const_params*r2*exp(-const_params*r2) - exp(-const_params*r2))/(r2*r2);

const_params2 = mueff*Power/(4.0*pi*(wPerf*cblood - kCond*mueff*mueff)*(wPerf*cblood - kCond*mueff*mueff));
D_1 = const_params2*(kCond*mueff*mueff*mueff*r1 + cblood*(2.0 - mueff*r1)*wPerf)*exp(-mueff*r1);
D_2 = const_params2*((kCond*mueff*mueff*mueff - cblood*mueff*wPerf)*exp(-mueff*r2) - mueff*(kCond*mueff*mueff*mueff*r2 + cblood*(2.0 - mueff*r2)*wPerf)*exp(-mueff*r2));
  //  printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n",mueff,Power,wPerf,cblood,kCond,mueff);
C_2 = (A_plus*D_2/B_plus - D_1)/(B_minus*A_plus/B_plus - A_minus);
C_1 = (D_1 - A_minus*C_2)/A_plus;

// Homogeneous solution
soln_h = C_1*exp(const_params*rVar)/rVar + C_2*exp(-const_params*rVar)/rVar;

// Particular solution
soln_p = const_params2*(kCond*mueff*mueff*mueff*rVar + cblood*(2.0 - mueff*rVar)*wPerf)*exp(-mueff*rVar)/rVar;

// General solution
sens = soln_h + soln_p;
*sensitivity = sens;

//printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n",mueff,Power,wPerf,soln_h,soln_p,soln_h+soln_p,sens,*sensitivity);
}
__device__
void pointSourceSensitivity(double rVar, double r1, double r2, double wPerf, double cblood, double kCond, double mueff, double u0, double ua, double Power, double *sensitivity)
{
double t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,
    t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,
    t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60;

  t2 = 1.0/sqrt(kCond);
  t3 = sqrt(cblood);
  t4 = sqrt(wPerf);
  t5 = sqrt(kCond);
  t6 = mueff*t5;
  t7 = t3*t4;
  t8 = t6+t7;
  t9 = 1.0/3.141592653589793;
  t10 = 1.0/sqrt(cblood);
  t11 = 1.0/sqrt(wPerf);
  t12 = exp(-rVar*t2*t8);
  t13 = 1.0/rVar;
  t14 = rVar*t2*t3*t4;
  t15 = t6-t7;
  t16 = t3*t4*2.0;
  t17 = mueff*mueff;
  t18 = exp(t14);
  t19 = cblood*wPerf;
  t20 = sqrt(t19);
  t21 = t6+t20;
  t22 = 1.0/t21;
  t23 = t6-t20;
  t24 = 1.0/t23;
  t25 = 1.0/(t8*t8);
  t26 = 1.0/(t15*t15);
  t27 = r1*t2*t3*t4*2.0;
  t28 = exp(t27);
  t29 = r2*t2*t3*t4*2.0;
  t30 = exp(t29);
  t31 = mueff*r1*t5;
  t32 = pow(cblood,3.0/2.0);
  t33 = pow(wPerf,3.0/2.0);
  t34 = mueff*r2*t5;
  t38 = r1*t3*t4;
  t35 = t34-t38;
  t36 = t2*t35;
  t37 = exp(t36);
  t39 = pow(kCond,3.0/2.0);
  t43 = r2*t3*t4;
  t40 = t31-t43;
  t41 = t2*t40;
  t42 = exp(t41);
  t44 = r2*r2;
  t45 = rVar*t3*t4;
  t46 = t19-kCond*t17;
  t47 = t5*t28;
  t48 = r2*t3*t4*t28;
  t49 = r2*t3*t4*t30;
  t50 = t47+t48+t49-t5*t30;
  t51 = 1.0/t50;
  t52 = cblood*cblood;
  t53 = wPerf*wPerf;
  t54 = t34+t38;
  t55 = t2*t54;
  t56 = exp(t55);
  t57 = t17*t17;
  t58 = t31+t43;
  t59 = t2*t58;
  t60 = exp(t59);
  *sensitivity = t13*t18*(mueff*t9*t10*t11*t12*t25*(t6+t16)*(1.0/8.0)-rVar*t9*t10*t11*t12*t17*t22*(1.0/8.0))-t5*t10*t11*t13*exp(-t14)*(mueff*t2*t9*t26*exp(-rVar*t2*t15)*(t6-t16)*(1.0/4.0)-rVar*t2*t9*t17*t18*t24*exp(-mueff*rVar)*(1.0/4.0))*(1.0/2.0)+mueff*t9*t10*t11*t13*t22*t24*t25*t26*t46*t51*exp(-t2*(t31+t34-t45))*(r2*t52*t53*t56*2.0+t5*t32*t33*t56*2.0-t5*t32*t33*t60*2.0-t20*t39*t44*t57*t60-mueff*r2*t5*t32*t33*t60*2.0+mueff*r1*t17*t20*t39*t56-cblood*mueff*r1*t5*t20*t56*wPerf+cblood*t5*t17*t20*t44*t60*wPerf-mueff*r1*r2*t20*t32*t33*t56+kCond*mueff*r1*r2*t3*t4*t17*t20*t56)*(1.0/4.0)+mueff*t9*t10*t11*t13*t22*t24*t25*t26*t46*t51*exp(-t2*(t31+t34+t45-r1*t3*t4*2.0-r2*t3*t4*2.0))*(r2*t37*t52*t53*2.0-t5*t32*t33*t37*2.0+t5*t32*t33*t42*2.0+t20*t39*t42*t44*t57-mueff*r1*t17*t20*t37*t39+mueff*r2*t5*t32*t33*t42*2.0+cblood*mueff*r1*t5*t20*t37*wPerf-cblood*t5*t17*t20*t42*t44*wPerf-mueff*r1*r2*t20*t32*t33*t37+kCond*mueff*r1*r2*t3*t4*t17*t20*t37)*(1.0/4.0);
//printf("%12.5e %12.5e %12.5e %12.5e %12.5e \n", rVar, t58,t59,t60, *sensitivity);
}
__device__
void DebugWrite(int idx,int idmat,double rad,double omega, double conduction, double mueff,double temp,double sens)
{
   printf("%d %d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",idx,idmat,rad,omega,conduction,mueff,temp,sens);
   //int j,k;

   //for (j=0;j<n;j++) {
   //   for (k=0;k<n+1;k++) {
   //      printf("%d %d %12.5e %12.5e ",k,j,a[k][j].real(),a[k][j].imag());
   //   }
   //   printf(" | %d  %12.5e %12.5e \n",j,x[j].real(),x[j].imag());
   //}
   //printf("\n");
}
/*
 * Device code
 */
__global__ 
void steadyStatePennesLaser(
         int const NTissue,
         const    int* MaterialID,
         const double* Perfusion,
         const double* ThermalConduction,
         const double* EffectiveAttenuation,
         double const innerRadius,
         double const outerRadius,
         int const NSource,
         double const Power,
         const double* SourceXloc,
         const double* SourceYloc,
         const double* SourceZloc,
         double const InitialTemperature,
         double const ArterialTemperature,
         double const SpecificHeatBlood,
	 double const SpacingX,
	 double const SpacingY,
	 double const SpacingZ,
         int const NpixelX,
         int const NpixelY,
         int const NpixelZ,
         double* d_TemperatureArray,
         double* d_SensitivityArray)
{

//     double SpacingX=0.00078;
    /*
      grid stride loop design pattern, 1-d grid
      http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
         - By using a loop, you can support any problem size even if it exceeds the largest grid size your CUDA device supports. Moreover, you can limit the number of blocks you use to tune performance.
    */
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; 
         idx < NpixelX * NpixelY * NpixelZ;
         idx += blockDim.x * gridDim.x) 
      {
        // compute indices
        int index = idx; // use dummy variable
        int kkk = index/(NpixelX*NpixelY); 
        index -= kkk*NpixelX*NpixelY; 
        
        int jjj = index/NpixelX; 
        index -= jjj*NpixelX; 
        
        int iii = index/1;

        /* get material parameters */
        int const idmaterial =  MaterialID[idx];
        double omega      = Perfusion[idmaterial];
        double conduction = ThermalConduction[idmaterial];
        double mueff      = EffectiveAttenuation[idmaterial];
	//printf("%d",mueff);
        // linear superpostion of temperature sources
        double temperature = 0.0;
        double sensitivity = 0.0;
        for (int lll=0;lll<NSource;lll++) 
          {
//           double radiusSQ = (iii * SpacingX + 0.13281 - SourceXloc[lll])*(iii * SpacingX + 0.13281 - SourceXloc[lll])
//                           + (jjj * SpacingY + 0.10547 - SourceYloc[lll])*(jjj * SpacingY + 0.10547 - SourceYloc[lll])
//                           + (kkk * SpacingZ + 0.06000 - SourceZloc[lll])*(kkk * SpacingZ + 0.06000- SourceZloc[lll]);

	   double radiusSQ=powf(iii*SpacingX-SourceXloc[lll],2)
			  +powf(jjj*SpacingY-SourceYloc[lll],2)
			  +powf(kkk*SpacingZ-SourceZloc[lll],2);//SourceXloc[0]*SourceXloc[0];
           double radius   = sqrt(radiusSQ);

           // call GF code 
	   double sourcetemperature;
           pointSource(radius, innerRadius, outerRadius, omega , SpecificHeatBlood, conduction , mueff, InitialTemperature, ArterialTemperature, Power , &sourcetemperature);

       double sourcesensitivity;
           pointSourceSensitivity(radius, innerRadius, outerRadius, omega, SpecificHeatBlood, conduction, mueff, InitialTemperature, ArterialTemperature, Power, &sourcesensitivity);

	   if (radius <= innerRadius && NSource ==1)
		{
                    sourcetemperature = InitialTemperature;
                    sourcesensitivity = 0;
printf("%12.5e %d Flag 1", radius, NSource);
		}
           if (radius <= innerRadius && NSource == 10)
		{
                    sourcetemperature = InitialTemperature+55;
                    sourcesensitivity = 0;
printf("%12.5e %d Flag 2", radius, NSource);
		}
           if (radius <= innerRadius && NSource > 1)
		{
                   sourcetemperature = InitialTemperature;
                   sourcesensitivity = 0;
//printf("%12.5e %d %12.5e %12.5e Flag 3 \n", radius, NSource, sourcetemperature, sourcesensitivity);
		}
           // DebugWrite(idx,idmaterial,radius,omega,conduction,mueff,sourcetemperature,sourcesensitivity);
           // superposition
	   if (idmaterial==0)
	 	{
		   temperature=0;
           sensitivity=0;
printf("%12.5e %d Flag 4", radius, NSource);
		}
	   else
		{
                   temperature = temperature + sourcetemperature/((double)NSource); 
                   sensitivity = sensitivity + sourcesensitivity/((double)NSource);
//printf("%12.5e %d Flag 5", radius, NSource);
		}	 
          }
        // store temperature in array
        d_TemperatureArray[idx] = temperature;
        d_SensitivityArray[idx] = sensitivity;
      }
}
