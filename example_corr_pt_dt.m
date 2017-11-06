%% This is a wrapper for invoking Oscar Miranda-Dominguez's corr_pt_dt.m for specific ROIs.
%% Written by Moosa Ahmed Oct 13, 2017
%%this code is built on Oscar Miranda-Dominguez original script templates.
%%
clear; clc;
%% Add the folders containing the function to your path
%addpath(genpath('/group_shares/fnl/bulk/code/utilities/corr_pt_dt'));
addpath(genpath('J:\code\internal\utilities\corr_pt_dt'));
%% folder path is to the folder above group_texts which contains paths to the processed directory of each subject sperated by group.
folder_path = 'I:\FAIR_HORAK\Projects\FOG_Oscar';
f = filesep;

%%define output ptseries folder name; i.e. name of parcellation 
%parc = 'R_PPN'
%%Select which groups paths you want to include
%MYlist = {'C_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','C_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','F_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','F_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','NF_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_210_skip_frames_5_TRseconds_2_50','NF_ALL_prisma_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_222_skip_frames_5_TRseconds_2_50','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
%MYlist = {'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
%MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_50_min_frames_165_skip_frames_5_TRseconds_2_00'};
MYlist = {'C_ALL_Trio_paths_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
%MYlist = {'F_ALL_Trio_paths-5-38_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00','NF_ALL_Trio_paths-1-29_MCMethod_power_2014_FD_only_FD_th_0_30_min_frames_159_skip_frames_5_TRseconds_2_00'};
parclist = {'R_6ma','R_6mp','L_6ma','L_6mp'}

for  parc = 1:length(parclist)
    for condition = 1:length(MYlist)
        txt_file = [folder_path f 'group_texts' f MYlist{condition} '_PC.txt'];
        FD_mask_name = [folder_path f 'Gui_envs\standard\Functional' f MYlist{condition} f 'frame_removal_mask.mat'];
        load(FD_mask_name)
        
        out_folder = [folder_path f 'Experiments' f 'Fconn' f parclist{parc} f MYlist{condition}];
        mkdir (out_folder)
        
        fid = fopen(txt_file);
        line1 = fgetl(fid);
        jj=0;
        while ischar(line1)
            try
                rp = [line1 f 'MNINonLinear/Results'];
                sub_id = rp((strfind(rp, 'EE_PD/')+6):(strfind(rp, '/20')-1));
                scan_date = rp((strfind(rp, '-SIEMENS')-8):(strfind(rp, '-SIEMENS')-1));
                
                dt_series_file = [rp f sub_id '_FNL_preproc_Atlas.dtseries.nii'];
                pt_series_file = [rp f sub_id '_FNL_preproc_SMA+STN+PPN_subcortical.ptseries.nii'];
                %pt_series_file = [rp f sub_id '_FNL_preproc_Atlas_FOX_ROI_combined.ptseries.nii'];
                %pt_series_file = [rp f sub_id '_FNL_preproc_Gordon.ptseries.nii'];
                %output_file = [out_folder f sub_id '_Gordon_' scan_date];
                output_file = [out_folder f sub_id '_' scan_date];
                if parclist{parc} == 'R_6ma'
                    ix = [1]';% indices of the ROIS in the dt_series to be used for correlation
                elseif parclist{parc} == 'R_6mp'
                    ix = [2]';
                elseif parclist{parc} == 'L_6ma'
                    ix = [3]';
                else
                    ix = [4]';
                end
                %ix = []
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