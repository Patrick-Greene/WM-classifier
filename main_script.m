% Script for running code from paper. Must first run extract_features.m to
% get features, which are saved in feature_data_path below.
clear; clc

test_subj_i = 1; % index of test subject to exclude from training data

% path to where feature data is saved (see extract_features_script.m)
feature_data_path = "/path to feature data";
max_ch = 16; % max number of channels on any shank
n_samples = 100; % number of samples to use in posterior predictive
alpha_cd = 1; % contact depth kernel width (alpha_2 in paper)


% load data and separate into train and test
all_subj_data = load(feature_data_path);
n_subj = length(all_subj_data);
train_subj_data = all_subj_data([1:test_subj_i-1, test_subj_1+1:n_subj]);
test_subj_data = all_subj_data(test_subj_i);

% get training set
[features_wm, features_gm, train_labels] = get_training_data(train_subj_data);

% get test set
[test_features, test_labels] = get_test_data(test_subj_data);

% Run the classifier
est_probs = wm_classifier(features_wm, features_gm, train_labels,...
    test_features, max_ch, alpha_cd, n_samples);


% Plot contact probabilities vs true labels and ROC curve
figure(1)
plot_marg_probs(est_probs, test_subj_data)

figure(2)
plot_ROC(est_probs, test_subj_data)












