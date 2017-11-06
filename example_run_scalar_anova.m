%% This is a wrapper for invoking Oscar Miranda-Dominguez's ciftify_scalar_anovan.m for specific ROIs.
%% Written by Moosa Ahmed July 11, 2017
%%this code is built on Oscar Miranda-Dominguez original script templates.
%%
clear; clc;
%% Set Paths
if ispc
    addpath(genpath('W:\code\development\utilities\fconn_anova'));
    addpath(genpath('W:\code\development\utilities\scalar_anova'));
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
folder_path = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments';
f = filesep;
MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};

%condition = 'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'

%%Rename for every new condition
for condition = 1:length(MYlist)
    Full_dt_series = [];
    dt_series_file = [folder_path f 'Fconn' f MYlist{condition} f '*.nii'];
    files = dir(dt_series_file);
    for file = files';
        dt_file = [folder_path f 'Fconn' f MYlist{condition} f file.name];
        cii = ciftiopen(dt_file,path_wb_c);
        newcii = cii;
        dt_series = double(newcii.cdata);
        Full_dt_series = horzcat(Full_dt_series,dt_series);
    end
    dt_name = [folder_path f 'Mat/Mat_files' f MYlist{condition} '.mat'];
    save(dt_name,'Full_dt_series');
end
%% Load different conditions .Mats 
C_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
F_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
NF_Trio_0_30 = struct2array(load('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Mat\Mat_files\NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00.mat'));
%% Horzcat and save
All_Trio_0_30 = horzcat(C_Trio_0_30, F_Trio_0_30, NF_Trio_0_30);
All_name = [folder_path f 'Mat/Mat_files' f 'All_Trio_0_30.mat']
save(All_name,'All_Trio_0_30');
%% Read dtseries
% cii=ciftiopen(dt_series_file,path_wb_c);
% newcii=cii;
% dt_series=double(newcii.cdata);
% 
% gy=64984;
% gy=3;
% np=120;
% n_times=4;
% ct=randn(gy,np);% data in format grayordinate x subject
% offset=0;
% for i=1:6
%     
%     foo=ct(:,(1:20)+offset)+10*i;
%     foo=foo.*repmat(1/n_times:1/n_times:1,gy,5);
%     
%     ct(:,(1:20)+offset)=foo;
%     offset=offset+20;
% end
% 
% % ct(:,1:end/2)=ct(:,1:end/2)+100;
% % ct(:,end/2+1:end)=ct(:,end/2+1:end)-100;
% % ct(:,1:40)=ct(:,1:40)*.05;
% % ct(:,81:120)=ct(:,81:120)*2;
% scalar_brainarea_subject=ct;
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
main_results = scalar_anovan(All_Trio_0_30,between_design,within_design,options)

toc
%% make ciftis
%path_wb_c='/usr/global/hcp_workbench/bin_linux64/wb_command'; %path to wb_command
%cifti_scalar_template='/group_shares/PSYCH/code/development/utilities/scalar_anova/sulc.32k_fs_LR.dscalar.nii';
cifti_scalar_template='I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Fconn\C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00/Az_307_20140627.dtseries.nii';
hosting_directory='I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\Analysis';

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

%%
filename='F_NF_kstest_p.dtseries.nii';
cii=ciftiopen(filename,path_wb_c);
newcii=cii;
foo=newcii.cdata;
foo=-log10(double(foo));

newcii.cdata=foo;


output_file=['minus_log10_' filename];
display(['Saving ' output_file]);
ciftisave(newcii,output_file,path_wb_c);

%% K-test
P_all=zeros(91282,1);
for i = 1:91282
    [H,P] = kstest2(F_Trio_0_30(i,:),NF_Trio_0_30(i,:));
    P_all(i)=P;
end

output_file = [hosting_directory f 'F_NF_kstest_p.dtseries.nii'];
cii=ciftiopen(cifti_scalar_template,path_wb_c);
newcii=cii;
newcii.cdata=P_all;
display(['Saving ' output_file]);
ciftisave(newcii,output_file,path_wb_c);
