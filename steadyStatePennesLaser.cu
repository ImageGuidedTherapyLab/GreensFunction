/*
 * Example Matlab cuda kernel interface.
 */


__device__
void pointSource(double r, double R1, double R2, double w, double k, double mueff, double u0, double ua, double P, double *temperature )
{

   double PI_Var = 3.141592653589793;
   *temperature = ua+(P*PI_Var*(mueff*mueff)*exp(-mueff*r)*(1.0/4.0))/(r*(w-k*(mueff*mueff)))-(exp(-R1*mueff-R2*mueff)*exp(r*sqrt(w/k))*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*R2*(mueff*mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-R1*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)))-(exp(R1*sqrt(w/k))*exp(R2*sqrt(w/k))*exp(-r*sqrt(w/k))*exp(-R1*mueff)*exp(-R2*mueff)*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)+P*PI_Var*R2*(mueff*mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)+R1*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)));
}
__device__
void DebugWrite(int idx,int idmat,double rad,double omega, double conduction, double mueff,double temp)
{
   printf("%d %d %12.5e %12.5e %12.5e %12.5e %12.5e\n",idx,idmat,rad,omega,conduction,mueff,temp);
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
         double const SpecificHeatBlood ,
         double const  SpacingX,
         double const  SpacingY,
         double const  SpacingZ,
         int const NpixelX,
         int const NpixelY,
         int const NpixelZ),
         double* d_TemperatureArray
{
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
        double omega      = SpecificHeatBlood * Perfusion[idmaterial];
        double conduction = ThermalConduction[idmaterial];
        double mueff      = EffectiveAttenuation[idmaterial];

        // linear superpostion of temperature sources
        double temperature = 0.0;
        for (int lll=0;lll<NSource;lll++) 
          {
           double radiusSQ = (iii * SpacingX - SourceXloc[lll])*(iii * SpacingX - SourceXloc[lll])
                           + (jjj * SpacingY - SourceYloc[lll])*(jjj * SpacingY - SourceYloc[lll])
                           + (kkk * SpacingZ - SourceZloc[lll])*(kkk * SpacingZ - SourceZloc[lll]);
           double radius   = sqrt(radiusSQ);
           // call GF code 
           double sourcetemperature;
           pointSource(radius, innerRadius, outerRadius, omega , conduction , mueff, InitialTemperature, ArterialTemperature, Power , &sourcetemperature);
           // DebugWrite(idx,idmaterial,radius,omega,conduction,mueff,sourcetemperature);
           // superposition
           temperature = temperature + sourcetemperature/((double)NSource); 
          }
        // store temperature in array
        d_TemperatureArray[idx] = temperature;
      }
}


