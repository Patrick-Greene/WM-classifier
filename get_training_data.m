% Organizes training data from training subset of all_subj_data
function [features_wm, features_gm, train_labels] = get_training_data(train_subj_data)
% make feature arrays; each is one long N x n_feats array concatenating the
% feature values for all patients and all contacts that are either gm or wm
% (excluding test_subj_i) 
[features_wm, features_gm] = make_feature_arrays(train_subj_data);

train_labels = {};
n_train_subj = length(train_subj_data);
for subj_i=1:n_train_subj
    if subj_i == excl_subj_i
        continue
    end
    n_shanks = length(train_subj_data{subj_i}.ch_type_list);
    for shank_i=1:n_shanks
        train_labels{end+1} = train_subj_data{subj_i}.ch_type_list{shank_i};
    end
end
end



function [features_wm, features_gm] = make_feature_arrays(train_subj_data)
% make feature arrays; each is N x n_feats array concatenating the
% feature values for all patients and all electrodes that are either gm or wm
% Set excl_subj_i > number of subjects to not exclude anyone
n_subj = length(train_subj_data);
features_wm = [];
features_gm = [];
for subj_i=1:n_subj
    n_shanks = length(train_subj_data{subj_i}.subj_data.ch_type_list);
    for shank_i=1:n_shanks
        wm_mask = logical(train_subj_data{subj_i}.subj_data.ch_type_list{shank_i});
        gm_mask = ~wm_mask;
        features_wm = [features_wm; train_subj_data{subj_i}.subj_data.ch_feat_list{shank_i}(wm_mask,:)];
        features_gm = [features_gm; train_subj_data{subj_i}.subj_data.ch_feat_list{shank_i}(gm_mask,:)];
    end
end

end
