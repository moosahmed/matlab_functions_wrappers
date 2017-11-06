%%build dscalars across subjects
%% Written by Moosa Ahmed Oct 10, 2017
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

%% upload dscalars

folder_path = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps';
f = filesep;

txt_file = [folder_path f 'made_dscalars_PC2.txt'];
%FD_mask_name = [folder_path f 'Gui_envs\standard\Functional' f MYlist{condition} f 'frame_removal_mask.mat'];
%load(FD_mask_name)

%out_folder = [folder_path f 'dscalars' f 'results'];
% mkdir (out_folder)

All_dscalars = []
fid = fopen(txt_file);
line1 = fgetl(fid);
jj=0;
while ischar(line1)
    try
        rp = [line1];
        %sub_id = rp((strfind(rp, 'EE_PD/')+6):(strfind(rp, '/20')-1));
        %scan_date = rp((strfind(rp, '-SIEMENS')-8):(strfind(rp, '-SIEMENS')-1));
        
        jj = jj+1;
        disp(jj);
        
        cii = ciftiopen(rp,path_wb_c);
        newcii = cii;
        dscalar = double(newcii.cdata);
        All_dscalars = horzcat(All_dscalars,dscalar);
        line1 = fgetl(fid);
    catch
    end
end
fclose(fid)

All_dscalars_full = All_dscalars(:,[1:17,35,36,38,44,46,51,58,64:67,87:110,18:34,37,39:43,45,47:50,52:57,59:63,111,68:86,112:120]);
%%
All_name = [folder_path f 'dscalars' f 'All_dscalars_full.mat']
save(All_name,'All_dscalars_full');
%%
All_dscalars_test = 1 - All_dscalars_full.^2;
All_dscalars_root = sqrt(All_dscalars_test);
%% Define the order of the indices
clear between_design
between_design(1).name='Dx';
between_design(1).subgroups(2).name='F';
between_design(1).subgroups(3).name='NF';
between_design(1).subgroups(1).name='C';
between_design(1).subgroups(1).ix=1:52;
between_design(1).subgroups(2).ix=53:92;
between_design(1).subgroups(3).ix=93:120;
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
main_results = scalar_anovan(atanh(All_dscalars_root),between_design,within_design,options)
toc
%% make ciftis
%path_wb_c='/usr/global/hcp_workbench/bin_linux64/wb_command'; %path to wb_command
%cifti_scalar_template='/group_shares/PSYCH/code/development/utilities/scalar_anova/sulc.32k_fs_LR.dscalar.nii';
cifti_scalar_template=['I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\dscalars/Az_307_20140627.dtseries.nii'];
hosting_directory=[folder_path f 'dscalars' f 'results' f 'full_group' f 'oneminus_rsquared'];
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
    [H,P] = kstest2(All_dscalars_test(i,1:52),All_dscalars_test(i,53:92));
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
    [H,P] = kstest2(All_dscalars_test(i,1:52),All_dscalars_test(i,93:120));
    P_all(i)=P;
end
output_file = [hosting_directory f 'F_C_kstest_p.dtseries.nii'];

cii=ciftiopen(cifti_scalar_template,path_wb_c);
newcii=cii;
newcii.cdata=P_all;
display(['Saving ' output_file]);
ciftisave(newcii,output_file,path_wb_c);
%% K-test NF C
P_all=zeros(91282,1);
for i = 1:91282
    [H,P] = kstest2(All_dscalars_test(i,53:92),All_dscalars_test(i,93:120));
    P_all(i)=P;
end
output_file = [hosting_directory f 'NF_C_kstest_p.dtseries.nii'];

cii=ciftiopen(cifti_scalar_template,path_wb_c);
newcii=cii;
newcii.cdata=P_all;
display(['Saving ' output_file]);
ciftisave(newcii,output_file,path_wb_c);
%% K-test Park C
P_all=zeros(91282,1);
for i = 1:91282
    [H,P] = kstest2(All_dscalars_test(i,1:92),All_dscalars_test(i,93:120));
    P_all(i)=P;
end
output_file = [hosting_directory f 'Park_C_kstest_p.dtseries.nii'];

cii=ciftiopen(cifti_scalar_template,path_wb_c);
newcii=cii;
newcii.cdata=P_all;
display(['Saving ' output_file]);
ciftisave(newcii,output_file,path_wb_c);
%%
%filename='F_C_kstest_p.dtseries.nii';
%filelist={'F_NF_kstest_p.dtseries.nii'};
%filename='F_NF_ttest_p.dtseries.nii';
%filelist = {'F_C_kstest_p.dtseries.nii','F_NF_kstest_p.dtseries.nii','NF_C_kstest_p.dtseries.nii'}
filelist = {'Park_C_kstest_p.dtseries.nii'};
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

%% Check group divergences
grouplist = {freezers,nonfreezers,controls}

for group = 1:length(grouplist)
    dscalar = tanh(mean(atanh(grouplist{group}),2));
    dscalar_test = 1-dscalar.^2;
    if group == 1
        name = 'freezers';
    elseif group == 2
        name = 'nonfreezers';
    else
        name = 'controls';
    end
    outfile = ['I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\dscalars\' name '_divergence_to_HCP.dtseries.nii'];
    %cii=ciftiopen(cifti_scalar_template,path_wb_c);
    newcii=cii;
    newcii.cdata=dscalar_test;
    display(['Saving ' outfile]);
    ciftisave(newcii,outfile,path_wb_c);
end




    