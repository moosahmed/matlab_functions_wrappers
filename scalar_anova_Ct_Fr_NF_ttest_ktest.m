%% This is a wrapper for invoking Oscar Miranda-Dominguez's ciftify_scalar_anovan.m for specific ROIs.
%% Written by Moosa Ahmed October 10th, 2017, also adding t-tests and kstests to this. 
%%this code is built on Oscar Miranda-Dominguez original script templates.
%%
clear; clc;
%% Set Paths
if ispc
    addpath(genpath('W:\code\development\utilities\fconn_anova'));
    addpath(genpath('W:\code\development\utilities\scalar_anova'));
    p{1}='W:/code/external/utilities/gifti-1.4';
    p{2}='W:/code/external/utilities/Matlab_CIFTI';
    for i=1:length(p)
        addpath(genpath(p{i}));
    end
else
    addpath(genpath('/group_shares/PSYCH/code/development/utilities/fconn_anova'));
    addpath(genpath('/group_shares/PSYCH/code/development/utilities/scalar_anova'));
    p{1}='/group_shares/PSYCH/code/external/utilities/gifti-1.4';
    p{2}='/group_shares/PSYCH/code/external/utilities/Matlab_CIFTI';
    for i=1:length(p)
        addpath(genpath(p{i}));
    end
end

%%Import Settings file
addpath('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat/support_files');
%addpath('/group_shares/fnl\code\internal\utilities\corr_pt_dt/support_files');
settings=settings_corr_pt_dt;%
%% Adding paths for this function
np=size(settings.path,2);
for i=1:np
    addpath(genpath(settings.path{i}));
end
path_wb_c=settings.path_wb_c; %path to wb_command
%% populate the data
folder_path = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\seed_vs_entire_brain';
f = filesep;
%parclist = {'R_6ma','R_6mp','L_6ma','L_6mp'}
%parclist = {'R_6ma','R_6mp'}
%parclist = {'L_6ma','L_6mp','L_PPN','R_PPN'}
%parclist = {'L_STN','L_PPN'}
parclist = {'Fox_14'}

%%
for parc = 1:length(parclist)
    %parc = 'R_STN'
    dt_folder = [folder_path f 'smoothed' f 'seed_to_brain_matrices_per_group' f parclist{parc}];
    mkdir (dt_folder)
    %%define output ptseries folder name; i.e. name of parcellation
    MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
    %MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
    %condition = 'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'
    %% Loop
    %%Rename for every new condition
    for condition = 1:length(MYlist)
        Full_dt_series = [];
        dt_series_file = [folder_path f 'smoothed' f 'seed_brain_correlation_ciftis_per_subj' f parclist{parc} f MYlist{condition} f '*.nii'];
        files = dir(dt_series_file);
        for file = files';
            dt_file = [folder_path f 'smoothed' f 'seed_brain_correlation_ciftis_per_subj' f parclist{parc} f MYlist{condition} f file.name];
            cii = ciftiopen(dt_file,path_wb_c);
            newcii = cii;
            dt_series = double(newcii.cdata);
            Full_dt_series = horzcat(Full_dt_series,dt_series);
        end
        dt_file = [dt_folder f MYlist{condition} '.mat'];
        save(dt_file,'Full_dt_series');
    end
    %% Load different conditions .Mats
    %C_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
    %F_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
    %NF_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
    C_Trio_0_30 = struct2array(load([folder_path f 'smoothed' f 'seed_to_brain_matrices_per_group' f parclist{parc} f 'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat']));
    F_Trio_0_30 = struct2array(load([folder_path f 'smoothed' f 'seed_to_brain_matrices_per_group' f parclist{parc} f 'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat']));
    NF_Trio_0_30 = struct2array(load([folder_path f 'smoothed' f 'seed_to_brain_matrices_per_group' f parclist{parc} f 'NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat']));
    %% Horzcat and save
    All_Trio_0_30 = horzcat(C_Trio_0_30, F_Trio_0_30, NF_Trio_0_30);
    All_name = [folder_path f 'smoothed' f 'seed_to_brain_matrices_per_group' f parclist{parc} f 'All_Trio_0_30.mat']
    save(All_name,'All_Trio_0_30');
    %%
    % Let's supose you have data from 3 groups
    %% Define the order of the indices
    clear between_design
    between_design(1).name='Dx';
    between_design(1).subgroups(1).name='C';
    between_design(1).subgroups(2).name='F';
    between_design(1).subgroups(3).name='NF';
    between_design(1).subgroups(1).ix=1:28;
    between_design(1).subgroups(2).ix=29:80;
    between_design(1).subgroups(3).ix=81:120;
    %between_design(2).name='Sex';
    %between_design(2).subgroups(1).name='Male';
    %between_design(2).subgroups(2).name='Female';
    %between_design(2).subgroups(1).ix=[1:20 41:60 81:100];
    %between_design(2).subgroups(2).ix=[21:40 61:80 101:120];
    %%
    within_design=[];
    options=[];
    %y=ct(1,:);
    %y=All_Trio_0_30(1,:);
    tic
    main_results = scalar_anovan(atanh(All_Trio_0_30),between_design,within_design,options)
    toc
    %% make ciftis
    %path_wb_c='/usr/global/hcp_workbench/bin_linux64/wb_command'; %path to wb_command
    %cifti_scalar_template='/group_shares/PSYCH/code/development/utilities/scalar_anova/sulc.32k_fs_LR.dscalar.nii';
    cifti_scalar_template=[folder_path f 'smoothed' f 'seed_brain_correlation_ciftis_per_subj' f parclist{parc} f 'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00/Az_307_20140627.dtseries.nii'];
    hosting_directory=[folder_path f 'smoothed' f 'Results_anova_k_t_test' f parclist{parc}];
    %%
    ciftify_scalar_anovan(main_results,path_wb_c,cifti_scalar_template,hosting_directory)
    %% Longitudinal
    % within_design(1).name='Time';
    % within_design(1).subgroups(1).name='m2';
    % within_design(1).subgroups(2).name='m4';
    % within_design(1).subgroups(3).name='m6';
    % within_design(1).subgroups(4).name='m8';
    % within_design(1).subgroups(1).ix=1:n_times:np;
    % within_design(1).subgroups(2).ix=2:n_times:np;
    % within_design(1).subgroups(3).ix=3:n_times:np;
    % within_design(1).subgroups(4).ix=4:n_times:np;
    %
    % main_results_long = scalar_anovan(scalar_brainarea_subject,between_design,within_design,options)
    %% make ciftis
    % path_wb_c='/usr/global/hcp_workbench/bin_linux64/wb_command'; %path to wb_command
    % cifti_scalar_template='/group_shares/PSYCH/code/development/utilities/scalar_anova/sulc.32k_fs_LR.dscalar.nii';
    % hosting_directory='/group_shares/FAIR_LAB/scratch/Oscar/setup_scalar_anova/test_ciftis_CT_long';
    %
    % ciftify_scalar_anovan(main_results_long,path_wb_c,cifti_scalar_template,hosting_directory)
    %% T-test
    [h,p] = ttest2(F_Trio_0_30',NF_Trio_0_30');
    output_file = [hosting_directory f 'F_NF_ttest_p.dtseries.nii'];
    
    cii=ciftiopen(cifti_scalar_template,path_wb_c);
    newcii=cii;
    newcii.cdata=p';
    display(['Saving ' output_file]);
    ciftisave(newcii,output_file,path_wb_c);
    %% K-test F NF
    P_all=zeros(91282,1);
    for i = 1:91282
        %[H,P] = kstest2(F_Trio_0_30(i,:),C_Trio_0_30(i,:));
        [H,P] = kstest2(F_Trio_0_30(i,:),NF_Trio_0_30(i,:));
        P_all(i)=P;
    end
    %output_file = [hosting_directory f 'F_C_kstest_p.dtseries.nii'];
    output_file = [hosting_directory f 'F_NF_kstest_p.dtseries.nii'];
    
    cii=ciftiopen(cifti_scalar_template,path_wb_c);
    newcii=cii;
    newcii.cdata=P_all;
    display(['Saving ' output_file]);
    ciftisave(newcii,output_file,path_wb_c);
    %% K-test F C
    P_all=zeros(91282,1);
    for i = 1:91282
        [H,P] = kstest2(F_Trio_0_30(i,:),C_Trio_0_30(i,:));
        P_all(i)=P;
    end
    output_file = [hosting_directory f 'F_C_kstest_p.dtseries.nii'];
    
    cii=ciftiopen(cifti_scalar_template,path_wb_c);
    newcii=cii;
    newcii.cdata=P_all;
    display(['Saving ' output_file]);
    ciftisave(newcii,output_file,path_wb_c);
    %%
    %filename='F_C_kstest_p.dtseries.nii';
    %filename='F_NF_kstest_p.dtseries.nii';
    %filename='F_NF_ttest_p.dtseries.nii';
    filelist = {'F_C_kstest_p.dtseries.nii','F_NF_kstest_p.dtseries.nii','F_NF_ttest_p.dtseries.nii'}
    for filename = 1:length(filelist)
        FULL_filename=[hosting_directory f filelist{filename}];
        cii=ciftiopen(FULL_filename,path_wb_c);
        newcii=cii;
        foo=newcii.cdata;
        foo=-log10(double(foo));
        
        newcii.cdata=foo;
        
        output_file=[hosting_directory f 'minus_log10_' filelist{filename}];
        display(['Saving ' output_file]);
        ciftisave(newcii,output_file,path_wb_c);
    end
    %% Controls Z-fissure
    C_atanh = atanh(C_Trio_0_30);
    M = mean(C_atanh,2);
    C_tanh = tanh(M);
    output_file = [hosting_directory f 'C_fiss.dtseries.nii'];
    cii=ciftiopen(cifti_scalar_template,path_wb_c);
    newcii=cii;
    newcii.cdata=C_tanh;
    display(['Saving ' output_file]);
    ciftisave(newcii,output_file,path_wb_c);
end