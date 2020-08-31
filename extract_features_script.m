% script to extract features from raw data and put into format required
% by classifier

% src_dir:  string giving location of patient data, i.e. "/home/my_data"
% save_dir: string giving location to save output
% file_names:  < n_patients, 1 > cell array of patient data file names (as strings)
% within src_dir
src_dir = "/path to data";
save_dir = "/path to save output";
file_names = {"subj_1", "subj_2"};


all_subj_data = {};
for subj_i=1:length(file_names)
    % load data into subj_data (do it like this bc load returns a
    % struct otherwise)
    subj_data = structfun(@(x) x, load(src_dir+"/"+file_names{subj_i}));
    all_subj_data{subj_i} = extract_features(subj_data);
end

save(save_dir+"/all_subj_data", 'all_subj_data')

