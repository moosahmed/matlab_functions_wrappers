%% This is a wrapper for invoking Robert Hermosillo and Mollie Marr's cifti_conn_matrixm and cifti_conn_pairwise_corr.m for specific ROIs.
%% Written by Moosa Ahmed September 20th, 2017
%%
clear all; clc; close all
%% run cifti_conn_matrix
if ispc
    tconc = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\15_smoothed_dt_paths_PC.conc';
    motion_conc = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\15_smoothed_motion_paths_PC.conc';
    conn_template = 'I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\scripts\from_mollie\Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
    addpath('I:\FAIR_HORAK\Projects\FOG_Oscar\Experiments\divergence_maps\scripts\from_mollie');
else
    tconc = '/group_shares/horaklab/bulk/FAIR_HORAK/Projects/FOG_Oscar/Experiments/divergence_maps/15_smoothed_dt_paths.conc';
    motion_conc = '/group_shares/horaklab/bulk/FAIR_HORAK/Projects/FOG_Oscar/Experiments/divergence_maps/15_smoothed_motion_paths.conc';
    conn_template = '/group_shares/horaklab/bulk/FAIR_HORAK/Projects/FOG_Oscar/Experiments/divergence_maps/scripts/from_mollie/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii';
    addpath('/group_shares/horaklab/bulk/FAIR_HORAK/Projects/FOG_Oscar/Experiments/divergence_maps/scripts/from_mollie');
end

series = 'd';
FD = 0.3;
TR = 1;
min = 2;
smoothing_kernal = 2.55;

%%
cifti_conn_matrix(tconc,strcat(series,'tseries'),motion_conc,FD,TR,min,smoothing_kernal)
subject_conn_conc = [tconc '_' series 'conn_of_' series 'tseries_' num2str(min) '_minutes_of_data_at_FD_' num2str(FD) '.conc'];
cifti_conn_pairwise_corr(conn_template,strcat(series,'conn'),subject_conn_conc)