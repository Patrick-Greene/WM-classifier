% script to extract features from raw data and put into format required
% by classifier
clear; clc

% Location of patient data, i.e. "/home/my_data"
src_dir = "/path to data";
% Location to save output
save_dir = "/path to save output";
% string array of patient data file names within src_dir
file_names = ["subj_1", "subj_2"]; 


for subj_i=1:length(file_names)
    "Extracting features for patient " + file_names(subj_i)
    % load data into subj_data (do it like this bc load returns a
    % struct otherwise)
    subj_data = structfun(@(x) x, load(src_dir+"/"+file_names{subj_i}));
    subj_data_features = extract_features(subj_data);
    
    save(save_dir+"/"+file_names(subj_i)+"_features", 'subj_data_features')
end

