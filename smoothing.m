%% this is a wb_command smoothing function. invoke using Run_smoothing.m
%% Written by Moosa Ahmed. Sept 13th, 2017

function output = smoothing(path_wb_c,dt_file,L_surf,R_surf,smoothing_kernal,direction)

out_file = [dt_file(1:length(dt_file)-13) '_SMOOTHED_' num2str(smoothing_kernal) '.dtseries.nii'];

if exist(out_file) ==0 %check to see if smoothed file exists yet.
    cmd = [path_wb_c ' -cifti-smoothing ' dt_file  ' ' smoothing_kernal ' ' smoothing_kernal ' ' direction ' ' out_file ' -left-surface ' L_surf ' -right-surface ' R_surf];
    system(cmd)
    %disp(cmd)
else %smoothed series already exists
    disp('Smoothed series already created for this subject')
end