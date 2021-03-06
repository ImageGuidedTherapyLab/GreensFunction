% Needs matrix in input

function [Temp,Sens] = GPU_kernel_calls_DM (power_log,spacing,scaling,mod_point,source,w_perf,k_cond,g_anisotropy,mu_eff,probe_u,robin_co,c_blood,choice)

%,g_anisotropy,mu_a,mu_s,probe_u,robin_co,c_blood,   matrix, power_log,source,w_perf,k_cond);

%% Embarrassing Parallel GPU Greens Function Linear Super Position
% clear all
% close all
format shortg


%% Simulate disjoint material/tissue types
% create npixel^3 image
% load rundata_rnd_6D_1e5_z
% load rundata_GL5_4D_low
if choice == 1            % mu
    Numruns = length(mu_eff);
    mu_eff=[zeros(Numruns,1) mu_eff'];
    k_cond=[0 k_cond];
    w_perf=[0 w_perf];
    
elseif choice ==2           % w_perf
    Numruns = length(w_perf);
    mu_eff = [0 mu_eff];
    if size(w_perf,1) > size(w_perf,2)
        w_perf = unique(w_perf);
        w_perf = w_perf';
        Numruns = length(w_perf);
    end
    w_perf=[zeros(Numruns,1) w_perf'];
    k_cond=[0 k_cond];
    
elseif choice ==3           % k_cond
    Numruns = length(k_cond);
    mu_eff = [0 mu_eff];
    k_cond=[zeros(Numruns,1) k_cond'];
    w_perf=[0 w_perf];
    
elseif choice ==4           % w_perf + mu
    Numruns = length(w_perf);
    mu_eff = [zeros(Numruns,1) mu_eff];
    w_perf = [zeros(Numruns,1) w_perf];
    k_cond=[0 k_cond];
    
elseif choice ==5           % w_perf + mu  RANDOM
    Numruns = length(w_perf);
    mu_eff = [zeros(Numruns,1) mu_eff'];
    k_cond=[0 k_cond];
    w_perf=[zeros(Numruns,1) w_perf'];
    
end



%para=[zeros(Numruns,1) para];
% npixel   = 100;
% materialID = int32(10*phantom3d('Modified Shepp-Logan',npixel));
% materialID(materialID == 3  ) = 1;
% materialID(materialID == 10 ) = 3;
% handle1 = figure(1)
% imagesc(materialID(:,:,npixel/2),[0 3])
% colorbar

%% Query the device
% GPU must be reset on out of bounds errors
reset(gpuDevice(1))
deviceInfo = gpuDevice(1);
numSMs = deviceInfo.MultiprocessorCount;


%% Setup Material Parameters
%ntissue=size(para,2);
ntissue=1;
% npixelx=domain.pointx;
% npixely=domain.pointy;
% npixelz=domain.pointz;
npixelx = mod_point.x*scaling.x;
npixely = mod_point.y*scaling.y;
npixelz = (mod_point.z * mod_point.z_subslice*scaling.z);
%npixelz_orig = npixelz;
npixelz = (npixelz + 1)/2;
%npixelz = (mod_point.z * mod_point.z_subslice*scaling.z) - (mod_point.z * mod_point.z_subslice*scaling.z -1 )/2 ;
% npixelx = mod_point.x;
% npixely = mod_point.y;
% npixelz = mod_point.z * mod_point.z;

%Tissue_region(:,:,:)=Tissue(domain.fov.x(1):domain.fov.x(2),domain.fov.y(1):domain.fov.y(2),domain.fov.z(1):domain.fov.z(2));
%materialID=int32(Tissue_region);
materialID = ones( npixelx, npixely, npixelz );
materialID = int32(materialID);
perfusion=w_perf;
conduction=k_cond;
% mueff=para(:,1:ntissue);
nsource=source.n;
% xloc=abs(laser.position(:,1)'-points.x(1));
% yloc=abs(laser.position(:,2)'-points.y(1));
% zloc=abs(laser.position(:,3)'-points.z(1));
xloc = source.laser(:,1);
yloc = source.laser(:,2);
zloc = source.laser(:,3);

u0 = probe_u - 37; % K
ua = 0; % K
u_artery=ua;
%c_blood=tissue.c_blood;
power=power_log;
R1=0.0015/2; % (m) R1 is the distance from the isotropic laser source point and the edge of the fiber
R2=1; % (m) R2 is the maximum edge of the domain;
spacingX=spacing.x;
spacingY=spacing.y;
spacingZ=spacing.z;

xloc = xloc + (spacingX * npixelx / 2);
yloc = yloc + (spacingY * npixely / 2);
%zloc = zloc + ( spacingZ * (npixelz / 2 - 1/2));

%% initialize data arrays
% initialize on host and perform ONE transfer from host to device
%npixelz=5;
h_temperature     = zeros(npixelx,npixely,npixelz);
d_temperature  = gpuArray( h_temperature  );
h_sensitivity = zeros(npixelx,npixely,npixelz);
d_sensitivity = gpuArray( h_sensitivity );

%% Compile and setup thread grid
% grid stride loop design pattern, 1-d grid
% http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
ssptx = parallel.gpu.CUDAKernel('steadyStatePennesLaserSensmod.ptx', 'steadyStatePennesLaserSensmod.cu');
%ssptx = parallel.gpu.CUDAKernel('SS_cuda1.ptx', 'SS_cuda1.cu');
threadsPerBlock = 256;
ssptx.ThreadBlockSize=[threadsPerBlock  1];
ssptx.GridSize=[numSMs*32               1];
% %% Run on GPU
% [d_temperature] = feval(ssptx,ntissue,materialID,perfusion,conduction, para(1,:), R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz, d_temperature);
% %%  transfer device to host
% h_temperature  = gather( d_temperature  );
% h_temperature=h_temperature+37;
% %%  plot temperature
% handle2 = figure(2)
% imagesc(h_temperature(:,:,1));
% colormap default
% colorbar

%tic
%[x0 y0 z0] = ndgrid(point
%For loop here
Temp = zeros(npixelx,npixely,Numruns);
Sens = zeros(npixelx,npixely,Numruns);

if choice ==1       % mu
    
    for ii = 1:Numruns
        ssptx = parallel.gpu.CUDAKernel('steadyStatePennesLaserSensmod.ptx', 'steadyStatePennesLaserSensmod.cu');
        %ssptx = parallel.gpu.CUDAKernel('SS_cuda1.ptx', 'SS_cuda1.cu');
        threadsPerBlock = 256;
        ssptx.ThreadBlockSize=[threadsPerBlock  1];
        ssptx.GridSize=[numSMs*32               1];
        
        if mod (ii,500) == 0
            toc
            fprintf('iter %d \n', ii);
        end
        %%  transfer device to host
        %         keyboard
        [d_temperature,d_sensitivity] = feval(ssptx,ntissue,materialID,perfusion,conduction, mu_eff(ii,:), R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz,d_temperature,d_sensitivity);
        tmp = gather( d_temperature );
        sen = gather( d_sensitivity );
        
        for kk = 1:npixelz-1
            tmp(:,:,npixelz+kk) = tmp(:,:,1+kk);
            sen(:,:,npixelz+kk) = sen(:,:,1+kk);
        end
        
        Temp(:,:,ii) = mean(tmp,3);
        Sens(:,:,ii) = mean(sen,3);
        if isnan(Temp)
            
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/troubleshoot
            save -v7.3 mu_nan.mat
            %   keyboard
        end
    end
    
    
    
elseif choice ==2   %w_perf
    
    for ii = 1:Numruns
        if mod (ii,500) == 0
            toc
            fprintf('iter %d \n', ii);
        end
        ssptx = parallel.gpu.CUDAKernel('steadyStatePennesLaserSensmod.ptx', 'steadyStatePennesLaserSensmod.cu');
        threadsPerBlock = 256;
        ssptx.ThreadBlockSize=[threadsPerBlock  1];
        ssptx.GridSize=[numSMs*32               1];
        %%  transfer device to host
        [d_temperature] = feval(ssptx,ntissue,materialID,perfusion(ii,:),conduction, mu_eff, R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz, d_temperature);
        tmp = gather( d_temperature );
        for kk = 1:npixelz-1
            tmp(:,:,npixelz+kk) = tmp(:,:,1+kk);
        end
        Temp(:,:,ii) = mean(tmp,3);
        if isnan(Temp)
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/troubleshoot
            save -v7.3 perf_nan.mat
            % keyboard
        end
    end
    
    
elseif choice ==3 %k_cond
    
    for ii = 1:Numruns
        if mod (ii,500) == 0
            toc
            fprintf('iter %d \n', ii);
        end
        ssptx = parallel.gpu.CUDAKernel('steadyStatePennesLaserSensmod.ptx', 'steadyStatePennesLaserSensmod.cu');
        threadsPerBlock = 256;
        ssptx.ThreadBlockSize=[threadsPerBlock  1];
        ssptx.GridSize=[numSMs*32               1];
        %%  transfer device to host
        [d_temperature] = feval(ssptx,ntissue,materialID,perfusion,conduction(ii,:), mu_eff, R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz, d_temperature);
        tmp = gather( d_temperature );
        for kk = 1:npixelz-1
            tmp(:,:,npixelz+kk) = tmp(:,:,1+kk);
        end
        Temp(:,:,ii) = mean(tmp,3);
        if isnan(Temp)
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/troubleshoot
            save -v7.3 cond_nan.mat
            % keyboard
        end
    end
    
elseif or(choice ==4 , choice ==5) %w_perf + mu
    
    for ii = 1:Numruns
        if mod (ii,500) == 0
            toc
            fprintf('iter %d \n', ii);
        end
        ssptx = parallel.gpu.CUDAKernel('steadyStatePennesLaserSensmod.ptx', 'steadyStatePennesLaserSensmod.cu');
        threadsPerBlock = 256;
        ssptx.ThreadBlockSize=[threadsPerBlock  1];
        ssptx.GridSize=[numSMs*32               1];
        %%  transfer device to host
        [d_temperature] = feval(ssptx,ntissue,materialID,perfusion(ii,:),conduction, mu_eff(ii,:), R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz, d_temperature);
        tmp = gather( d_temperature );
        for kk = 1:npixelz-1
            tmp(:,:,npixelz+kk) = tmp(:,:,1+kk);
        end
        Temp(:,:,ii) = mean(tmp,3);
        if isnan(Temp)
            cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/troubleshoot
            save -v7.3 cond_nan.mat
            % keyboard
        end
    end
    
end
toc
%keyboard
end


%
% keyboard
% %%  global search and plot exhaustive search
% tic
% Temp=zeros(numel(points.x),numel(points.y),5,Numruns);
% [x0 y0 z0]=ndgrid(points.x,points.y,points.z);
% [x y z]=ndgrid(points.x,points.y,[0.06 0.056 0.052 0.048 0.044]);
% % Tmap=34.6*ones(256,256,5);
% tic
% V=readVTK('dog1dtmap.0077.vtk');
% Tm=V(139:165,100:125,:);
% Tmeas(:,:,1)=Tm(:,:,1)';
% Tmeas(:,:,2)=Tm(:,:,2)';
% % Tmeas(:,:,3)=Tm(:,:,3)';
% Tmeas(:,:,4)=Tm(:,:,4)';
% Tmeas(:,:,5)=Tm(:,:,5)';
% Err=inf(Numruns,1);
% for iii = 1:Numruns
%     if mod(iii,100 )==0
%      toc
%     disp(sprintf('iter %d',iii));
%     end
%     %  h_temperature=zeros(numel(points.x),numel(points.y),numel(points.z));
%     %  d_temperature=gpuArray( h_temperature  );
%      [d_temperature] = feval(ssptx,ntissue,materialID,perfusion,conduction, para(iii,:), R1, R2, nsource, power ,xloc,yloc,zloc, u0 ,u_artery , c_blood, spacingX,spacingY,spacingZ,npixelx,npixely,npixelz, d_temperature);
%      tmp  = gather( d_temperature  );
%      tmp=tmp+34.6;
%     %  imagesc(tmp(:,:,10))
%      Ttmp=interpn(x0,y0,z0,double(tmp),x,y,z);
%      Temp(:,:,:,iii)=Ttmp;
%       Ttmp1=Ttmp(11:36,10:36,:);
%      Ttmp1(:,:,3)=0;
% %      norm(Ttmp(:)-Tmeas(:))
%      Err(iii)=norm(Ttmp1(:)-Tmeas(:));
% %      min(Err)
%     %  Tmap(domain.fov.x(1):domain.fov.x(2),domain.fov.y(1):domain.fov.y(2),:)=Ttmp;
% end
% toc
