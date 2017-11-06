%% This script creates smoothed dconns, parcellates them, then invokes Oscar Miranda-Dominguez's corr_pt_dt.m function
%% Written by Moosa Ahmed September 26th, 2017
%%
clear; clc;
%% Add the folders containing the function to your path
%addpath(genpath('/group_shares/fnl/bulk/code/utilities/corr_pt_dt'));
addpath(genpath('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\dorsal_medial_cerebellum_vs_cortex'));
%% Add the folders containing the function to your path
%addpath(genpath('/group_shares/fnl/bulk/code/utilities/corr_pt_dt'));
addpath(genpath('J:\code\internal\utilities\corr_pt_dt'));
%% Set Paths
if ispc
    p{1}='W:/code/external/utilities/gifti-1.4';
    p{2}='W:/code/external/utilities/Matlab_CIFTI';
    for i=1:length(p)
        addpath(genpath(p{i}));
    end
else
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
%% folder path is to the folder above group_texts which contains paths to the processed directory of each subject sperated by group.
folder_path = 'I:\FAIR_HORAK\Projects\FOG_Oscar';
parcel_path = 'I:/FAIR_HORAK/parcellations';
f = filesep;
smoothing_kernal = '2.55';

%%define output ptseries folder name; i.e. name of parcellation 
%parc = 'R_PPN'

%%Select which groups paths you want to include
%MYlist = {'C_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','C_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','F_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','F_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','NF_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','NF_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
%MYlist = {'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
%MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
%MYlist = {'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};

for condition = 1:length(MYlist)
    txt_file = [folder_path f 'group_texts' f MYlist{condition} '_PC.txt'];
    
    fid = fopen(txt_file);
    line1 = fgetl(fid);
    while ischar(line1)
        try
            rp = [line1 f 'MNINonLinear/Results'];
            sub_id = rp((strfind(rp, 'EE_PD/')+6):(strfind(rp, '/20')-1));
            scan_date = rp((strfind(rp, '-SIEMENS')-8):(strfind(rp, '-SIEMENS')-1));
            
            dt_file = [rp f sub_id '_FNL_preproc_Atlas.dtseries.nii'];
            L_surf = [line1 f 'MNINonLinear/fsaverage_LR32k' f sub_id '.L.midthickness.32k_fs_LR.surf.gii'];
            R_surf = [line1 f 'MNINonLinear/fsaverage_LR32k' f sub_id '.R.midthickness.32k_fs_LR.surf.gii'];
                
            disp(dt_file);
            %disp(L_surf);
            %disp(R_surf);
            
            smoothing(path_wb_c,dt_file,L_surf,R_surf,smoothing_kernal,'COLUMN')
            
            smoothed_file = [dt_file(1:length(dt_file)-13) '_SMOOTHED_' num2str(smoothing_kernal) '.dtseries.nii'];
            
            parcelate = [path_wb_c ' -cifti-parcellate ' smoothed_file ' ' parcel_path f 'FOX_ROI_combined.dlabel.nii' ' COLUMN ' rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_FOX_ROI_combined.ptseries.nii'];
            parcelate2 = [path_wb_c ' -cifti-parcellate ' smoothed_file ' ' parcel_path f 'SMA+STN+PPN.32k_fs_LR.dlabel.nii' ' COLUMN ' rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_SMA+STN+PPN.32k_fs_LR.ptseries.nii'];
            
            system(parcelate)
            system(parcelate2)
            
            line1 = fgetl(fid);
        catch
        end
    end
    fclose(fid)
end

% %% parcelate
% for condition = 1:length(MYlist)
%     txt_file = [folder_path f 'group_texts' f MYlist{condition} '_PC.txt'];
%     
%     fid = fopen(txt_file);
%     line1 = fgetl(fid);
%     while ischar(line1)
%         try
%             rp = [line1 f 'MNINonLinear/Results'];
%             sub_id = rp((strfind(rp, 'EE_PD/')+6):(strfind(rp, '/20')-1));
%             scan_date = rp((strfind(rp, '-SIEMENS')-8):(strfind(rp, '-SIEMENS')-1));
%             
%             dt_file = [rp f sub_id '_FNL_preproc_Atlas.dtseries.nii'];
%                 
%             disp(dt_file);
%             %disp(L_surf);
%             %disp(R_surf);
%             
%             smoothed_file = [dt_file(1:length(dt_file)-13) '_SMOOTHED_' num2str(smoothing_kernal) '.dtseries.nii'];
%             
%             parcelate = [path_wb_c ' -cifti-parcellate ' smoothed_file ' ' parcel_path f 'FOX_ROI_combined.dlabel.nii' ' COLUMN ' rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_FOX_ROI_combined.ptseries.nii'];
%             parcelate2 = [path_wb_c ' -cifti-parcellate ' smoothed_file ' ' parcel_path f 'SMA+STN+PPN.32k_fs_LR.dlabel.nii' ' COLUMN ' rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_SMA+STN+PPN.32k_fs_LR.ptseries.nii'];
%             
%             system(parcelate)
%             system(parcelate2)
%             
%             line1 = fgetl(fid);
%         catch
%         end
%     end
%     fclose(fid)
% end
%% Corr_pt_dt
%parclist = {'R_6ma','R_6mp','L_6ma','L_6mp','L_PPN','R_PPN','L_STN','R_STN'}
parclist = {'Fox_14'}

for  parc = 1:length(parclist)
    for condition = 1:length(MYlist)
        txt_file = [folder_path f 'group_texts' f MYlist{condition} '_PC.txt'];
        FD_mask_name = [folder_path f 'Gui_envs\standard\Functional' f MYlist{condition} f 'frame_removal_mask.mat'];
        load(FD_mask_name)
        
        out_folder = [folder_path f 'Experiments' f 'Fconn' f 'smoothed' f parclist{parc} f MYlist{condition}];
        mkdir (out_folder)
        
        fid = fopen(txt_file);
        line1 = fgetl(fid);
        jj=0;
        while ischar(line1)
            try
                rp = [line1 f 'MNINonLinear/Results'];
                sub_id = rp((strfind(rp, 'EE_PD/')+6):(strfind(rp, '/20')-1));
                scan_date = rp((strfind(rp, '-SIEMENS')-8):(strfind(rp, '-SIEMENS')-1));
                
                dt_series_file = [rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55.dtseries.nii'];
                %pt_series_file = [rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_SMA+STN+PPN.32k_fs_LR.ptseries.nii'];
                pt_series_file = [rp f sub_id '_FNL_preproc_Atlas_SMOOTHED_2.55_FOX_ROI_combined.ptseries.nii'];
                %pt_series_file = [rp f sub_id '_FNL_preproc_Atlas_FOX_ROI_combined.ptseries.nii'];
                %pt_series_file = [rp f sub_id '_FNL_preproc_Gordon.ptseries.nii'];
                %output_file = [out_folder f sub_id '_Gordon_' scan_date];
                output_file = [out_folder f sub_id '_' scan_date];
%                 if parclist{parc} == 'R_6mp'
%                     ix = [2]';% indices of the ROIS in the dt_series to be used for correlation
%                 elseif parclist{parc} == 'L_6ma'
%                     ix = [3]';
%                 elseif parclist{parc} == 'L_6mp'
%                     ix = [4]';
%                 elseif parclist{parc} == 'L_PPN'
%                     ix = [5]';
%                 elseif parclist{parc} == 'R_PPN'
%                     ix = [6]';
%                 elseif parclist{parc} == 'L_STN'
%                     ix = [7]';
%                 elseif parclist{parc} == 'R_STN'
%                     ix = [8]';
%                 else
%                     ix = [1]';
%                 end
                ix = [1]';
                jj = jj+1;
                FD_mask = mask{jj};
                
                disp(dt_series_file);
                disp(pt_series_file);
                disp(jj);
                disp(length(FD_mask));
                
                corr_pt_dt(pt_series_file,dt_series_file,ix,output_file,FD_mask)
                line1 = fgetl(fid);
            catch
            end
        end
        fclose(fid)
    end
end