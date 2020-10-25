% Script for running code from paper. Must first run extract_features.m to
% get features, which are saved in feature_data_path below.


clear; clc

% path to where feature data is saved (see extract_features_script.m)
feature_data_dir = "/home/pg/my_programs/sarma_lab/bids_layout_data/wm_gm_features2";
file_names = "subj_" + [1:29]; % file names of feature data
% created from extract_features_script.m
n_folds = 4; % number of cross val folds
max_ch = 16; % max number of channels on any shank
n_samples = 100; % number of samples to use in posterior predictive
alpha_cd = 1; % contact depth kernel width (alpha_2 in paper)

n_subj = length(file_names);
partitions = cvpartition(n_subj,'KFold', n_folds); % 4 fold cross val

AUC = [];
for fold_i=1:n_folds
    test_subj_indices = find(test(partitions, fold_i));
    train_subj_indices = find(training(partitions, fold_i));

    % load data for train and test sets
    train_subj_data = {};
    for i=train_subj_indices'
        % have to do this stupid thing bc matlab load returns a struct
        temp=load(feature_data_dir+"/"+file_names(i));
        field = fieldnames(temp);
        train_subj_data{end+1} = temp.(field{1});
    end

    test_subj_data = {};
    for i=test_subj_indices'
        temp=load(feature_data_dir+"/"+file_names(i));
        field = fieldnames(temp);
        test_subj_data{end+1} = temp.(field{1});
    end


    % format training set
    [features_wm, features_gm, train_labels] = get_training_data(train_subj_data);
    
    for i=1:length(test_subj_data)
        % format test set
        [test_features, test_labels] = get_test_data(test_subj_data{i});

        % Run the classifier
        est_probs = wm_classifier(features_wm, features_gm, train_labels,...
            test_features, max_ch, alpha_cd, n_samples);
        
        AUC(end+1) = plot_ROC(est_probs, test_subj_data{i});
        AUC(end+1)
    end
end

mean(AUC)
std(AUC)













