% This script finds the best mu_eff for the different studies.
function init_and_save (choice);
%choice = 1; % 1 = mu; 2 = perf; 3 = cond;
tic

%cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
% data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
% 
% %opttype = 'bestfit50';
% datasummary = dlmread(data_filename,',',1,0);
% datasummary(any(isnan(datasummary), 2), 7) = 1;
% num_studies = size(datasummary,1);
% for ii = 1:num_studies
%     
%     Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
%     Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
%     
% end
% clear ii

%[~, ~, max_phys_sz]=find_spatial_dim(1);
FOV = [0.04 0.04]; % FOV in meters
min_space = 0.00085;
num_pix = FOV./min_space;
max_phys_sz = [FOV(1), min_space, min_space, num_pix(1); FOV(2), min_space, min_space, num_pix(2)];

%% choice = 1 is mu; 2 is perf; 3 is cond;
% input_path = cell(1,2);
% input_path{1,1} = Study_paths{1,1};
% input_path{1,2} = Study_paths{1,2};
[ all_opt_fig, no_pwr_fig,sim_dim, summary ] = set_mu_values ( max_phys_sz, choice );

cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries

if choice == 1
    
    %save ('all_opt_mu.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    save ('all_opt_mu_small.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    
elseif choice == 2
    
    save ('all_opt_perf.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    
elseif choice == 3
    
    save ('all_opt_cond.mat','all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    
elseif choice ==4
    
    %save ('all_opt_perf_mu_400_short.mat', 'all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    save ('opt_perf_mu_400_long.mat', 'all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
    
elseif choice ==5
    
    save ('rand_opt_perf_mu222222.mat', 'all_opt_fig','no_pwr_fig','sim_dim','summary','-v7.3');
end
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
end