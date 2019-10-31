data = readtable('ADHD_ICM_random200.xlsx');

structures = ["x_Amyg_R_","x_Fimbria_R_","x_Hippo_R_","x_Mammillary_R_","x_Amyg_L_","x_Fimbria_L_","x_Hippo_L_","x_Mammillary_L_"];
%%

% group controls and patients
% 108 ADHD rows 
sorted_data = sortrows(data,4);

% extract structures of interest
for i = 1:length(structures)
    our_struc(:,i) = sorted_data(:,char(structures(i)));
end

% split controls and patients stuctures
patient_struc = our_struc(1:108,:);
control_struc = our_struc(109:end,:);

% convert table to array
patient_array = table2array(patient_struc);
control_array = table2array(control_struc);
all_array = table2array(our_struc);

% find mean and stddev of each structure for controls
vol_mean = mean(control_array,1);
vol_stddev = std(control_array);

% structure loop
for i = 1:size(all_array,2)
    % each subject loop
	for j = 1:size(all_array,1)
        zscore(j,i) = (all_array(j,i)-vol_mean(i))/vol_stddev(i);
    end
end

a = heatmap(zscore(1:30,:));

% % convert zscore array to table
% zscore_table = array2table(zscore);
% zscore_table.Properties.VariableNames = {'Amyg_R','Fimbria_R','Hippo_R','Mammillary_R','Amyg_L','Fimbria_L','Hippo_L','Mammillary_L'};

