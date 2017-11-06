%% Makes 3D Structural Matrices ROIxROIxSubject
%% Written By Moosa Ahmed Oct 1st,2017
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
%% Upload Struc Matrices
folder_path = 'Y:\Projects\Allen-HumanAdult-OHSU\struc_from_alan\N80\N80';
mat_files_folder_path = 'Y:\Projects\Allen-HumanAdult-OHSU\struc_from_alan\Mat_files\N31';
f = filesep;
%%
%MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
%MYlist = {'107321','107422'};
connlist = {'Conn3','Conn1'};
% MYlist = {'105216',
%     '107321',
%     '107422',
%     '110411',
%     '139233',
%     '146432',
%     '158035',
%     '159340',
%     '164030',
%     '173435',
%     '173940',
%     '185442',
%     '194645',
%     '198855',
%     '201111',
%     '208327',
%     '214221',
%     '268850',
%     '293748',
%     '386250',
%     '510326',
%     '567961',
%     '573249',
%     '627549',
%     '654754',
%     '695768',
%     '789373',
%     '833148',
%     '837560',
%     '877269',
%     '912447'};
%condition = 'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'
%% Loop
%%Rename for every new condition
for conn = 1:length(connlist);
    Full_struc = [];
    struc_file = [folder_path f connlist{conn} f '*.pconn.nii'];
    files = dir(struc_file);
    for file = files';
        S_file = [folder_path f connlist{conn} f file.name];
        cii = ciftiopen(S_file,path_wb_c);
        newcii = cii;
        struc = double(newcii.cdata);
        Full_struc = cat(3,Full_struc,struc);
    end
    struc_name = [mat_files_folder_path f connlist{conn} '.mat'];
    save(struc_name,'Full_struc');
end

% for condition = 1:length(MYlist)
%         Full_dt_series = [];
%         dt_series_file = [folder_path f 'smoothed' f 'seed_brain_correlation_ciftis_per_subj' f parclist{parc} f MYlist{condition} f '*.nii'];
%         files = dir(dt_series_file);
%         for file = files';
%             dt_file = [folder_path f 'smoothed' f 'seed_brain_correlation_ciftis_per_subj' f parclist{parc} f MYlist{condition} f file.name];
%             cii = ciftiopen(dt_file,path_wb_c);
%             newcii = cii;
%             dt_series = double(newcii.cdata);
%             Full_dt_series = horzcat(Full_dt_series,dt_series);
%         end
%         dt_file = [dt_folder f MYlist{condition} '.mat'];
%         save(dt_file,'Full_dt_series');
%     end