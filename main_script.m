% Script for running code from paper. Must first run extract_features.m to
% get features, which are saved in feature_data_path below.
clear; clc

test_subj_i = 1; % index of test subject to exclude from training data

% path to where feature data is saved (see extract_features_script.m)
feature_data_dir = "/path to feature data";
file_names = ["subj_1", "subj_2"]+"_features"; % file names of feature data
% created from extract_features_script.m

max_ch = 16; % max number of channels on any shank
n_samples = 100; % number of samples to use in posterior predictive
alpha_cd = 1; % contact depth kernel width (alpha_2 in paper)


% load data for train and test sets
n_subj = length(file_names);
train_subj_data = {};
for i=[1:test_subj_i-1, test_subj_i+1:n_subj]
    % have to do this stupid thing bc matlab load returns a struct
    temp=load(feature_data_dir +"/"+ file_names(i));
    field = fieldnames(temp);
    train_subj_data{end+1} = temp.(field{1});
end

test_subj_data = {};
temp=load(feature_data_dir +"/"+ file_names(test_subj_i));
field = fieldnames(temp);
test_subj_data{1} = temp.(field{1});


% format training set
[features_wm, features_gm, train_labels] = get_training_data(train_subj_data);

% format test set
[test_features, test_labels] = get_test_data(test_subj_data{1});

% Run the classifier
est_probs = wm_classifier(features_wm, features_gm, train_labels,...
    test_features, max_ch, alpha_cd, n_samples);


% Plot contact probabilities vs true labels and ROC curve
figure(1)
plot_marg_probs(est_probs, test_subj_data{1})

figure(2)
plot_ROC(est_probs, test_subj_data{1})












