% Organizes test data from test set subset of all_subj_data
function [test_features, test_labels] = get_test_data(test_subj_data)
% test_features{shank_i} is the feature matrix (n_elecs x n_features) for shank_i
% test_labels{shank_i} are the true labels (1=wm, 0=gm) corresponding to shank_elecs{shank_i}

test_features = {};
test_labels = {};
n_shanks = length(test_subj_data{1}.ch_type_list);
for shank_i=1:n_shanks
    test_features{shank_i} = test_subj_data{1}.ch_feat_list{shank_i};
    test_labels{shank_i} = test_subj_data{1}.ch_type_list{shank_i};
end
end
