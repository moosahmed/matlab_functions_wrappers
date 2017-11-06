%% this code populates gene expression values into the subject's surface
%% Code written by Oscar Miranda-Dominguez and Moosa Ahmed May 9th, 2017

%% Add code to your working session

if ispc
    p{1}='P:\code/external/utilities/gifti-1.6';
    p{2}='P:\code/external/utilities/Matlab_CIFTI';
    cd M:\scratch\Oscar\ciftis\coordinates
else
    p{1}='/group_shares/PSYCH/code/external/utilities/gifti-1.6';
    p{2}='/group_shares/PSYCH/code/external/utilities/Matlab_CIFTI';
    cd /group_shares/FAIR_LAB/scratch/Oscar/ciftis/coordinates/
end
for i=1:length(p)
    addpath(genpath(p{i}));
end
%% Edit on March 16 to mess with func.gii
% 1. Define the path to the func gii file that will be used as template.
% Make sure this file is within ~/MNINonLinear/fsaverage_LR32k for the
% particpant
%
% There is one foer each hemisphere
if ispc
    %func_file='T:\Projects\Allen-HumanAdult-OHSU\scratch_TEMP_NO_BACKUP\HCP_Processing\H2002\HCP_release_20160610_NoT2\H2002\MNINonLinear\fsaverage_LR32k\H2002.L.MyelinMap.32k_fs_LR.func.gii';
else
	func_file='/group_shares/FAIR_HCP/Projects/Allen-HumanAdult-OHSU/scratch_TEMP_NO_BACKUP/New_HCP_Processing/H1015/HCP_release_20161027/H1015/MNINonLinear/fsaverage_LR32k/H1015.L.MyelinMap.32k_fs_LR.func.gii';
end

%% 2. Read the content of the file into g
g = gifti(func_file);
ng=size(g.cdata,1); 
%% 
ix_file='/group_shares/FAIR_LAB/scratch/Oscar/ciftis/coordinates/H1015/H1015_L_vertex_list_out.txt';% This file lists the grayordinates associated to each probe. It was made by wb_command
%ix_file='M:/scratch/Oscar/ciftis/coordinates/H2002/MNI_H2002_L_vertex_list_out.txt';% This file lists the grayordinates associated to each probe. It was made by wb_command
fileID = fopen(ix_file);
C = textscan(fileID,'%d');
fclose(fileID);

%ix=C{1};
ix=C{1}+1;%The standard index space is zero indexed. We are adding and offset of one because the vectors in matlab are base-1

%%
for ii = 0:16
	j = num2str(ii);
	filename=['/group_shares/FAIR_HCP/Projects/Allen-HumanAdult-OHSU/scratch_TEMP_NO_BACKUP/normalized_microarray_donor15496/for_dt/L/' j '.csv'];
	%'/group_shares/FAIR_LAB/scratch/Oscar/ciftis/coordinates/H2002/trimd_MicroarrayExpression_cortex.csv'
	%gen_data = import_gene_data_file(filename);
	gen_data = csvread(filename);
    [r c]=size(gen_data);

	temp=zeros(ng,r);
	temp(ix,:)=gen_data';

	gg=g;
	gg.cdata=temp;
	save(gg,['/group_shares/FAIR_HCP/Projects/Allen-HumanAdult-OHSU/dtseries_Files/Funcgii/H1015/L/' j '-H1015.L.Gene.32k_fs_LR.func.gii']);
end
% %% Define the location of wb_command
%path_wb_c='/usr/global/hcp_workbench/bin_linux64/wb_command'; %path to wb_command
%% Read a dtseries to be used as template that will be replace using gene data
% You can read label files (parcels like gordon), dense time courses,
% parcellated time courses, etc

% Define the file you would like to read
% file_path='/group_shares/FAIR_LAB/scratch/Oscar/ciftis';

% Pick one of the following files to see the content
% filename='10050-2_FNL_preproc_Atlas.dtseries.nii';
%filename='10050-2_FNL_preproc_Gordon.ptseries.nii';
% filename='10050-2_FNL_preproc_Gordon_subcortical.ptseries.nii';
% filename='Gordon_333_Parcel+FS_Anatomical_Parcel.dlabel.nii';
% file=[file_path '/' filename];

%cii=ciftiopen(file,path_wb_c);
%newcii=cii;
%X=newcii.cdata;
%whos X
%% Get the number of grayordinates to preallocate memory later
%[r c]=size(X);

%% Read gene expresion
%for ii = 0:16
    %j = num2str(ii)
    %filename=['/group_shares/FAIR_HCP/Projects/Allen-HumanAdult-OHSU/scratch_TEMP_NO_BACKUP/normalized_microarray_donor10021/for_dt/L/' j '.csv'];
    %'/group_shares/FAIR_LAB/scratch/Oscar/ciftis/coordinates/H2002/trimd_MicroarrayExpression_cortex.csv'
    %gen_data = import_gene_data_file(filename);
    
    %[genes,probes]=size(gen_data);
    %Y=zeros(r,genes);
    
    %ix_file='/group_shares/FAIR_LAB/scratch/Oscar/ciftis/coordinates/H2002/MNI_H2002_L_vertex_list_out.txt';% This file lists the grayordinates associated to each probe. It was made by wb_command
    %fileID = fopen(ix_file);
    %C = textscan(fileID,'%d');
    %fclose(fileID);
    
    %ix=C{1};
    %ix=C{1}+32492;
    %Y(ix,:)=gen_data';% NOte the transpose since the original csv file has in rows the number of genes, and the probes are coded as columns
    %newcii.cdata=Y;
    %% Populate gene expresion as cifti
    %output_file=['/group_shares/FAIR_HCP/Projects/Allen-HumanAdult-OHSU/dtseries_Files/H2002/L/' j 'MNI.dtseries.nii']; % defining the name of your output file
    %ciftisave(newcii,[output_file],path_wb_c); % Making your cifti
    %system(['gzip ' output_file]);
%end
