function AUC = plot_ROC(est_probs, test_subj_data)

thresh_to_use = [1:-0.01:0];
test_labels = test_subj_data.ch_type_list;
n_thresh = length(thresh_to_use);
false_pos = nan(1,n_thresh);
true_pos = nan(1,n_thresh);
for thresh_i=1:n_thresh
    [false_pos(thresh_i),true_pos(thresh_i)] = calc_fp_tp(thresh_to_use(thresh_i), est_probs, test_labels);
end
[fp5, tp5] = calc_fp_tp(0.5, est_probs, test_labels); % fp,tp rate at thresh=0.5


% Get area under ROC curve
AUC = calc_AUC(false_pos, true_pos);

hold on
plot(false_pos, true_pos,'k')
plot(thresh_to_use, thresh_to_use, 'k--')

scatter(fp5,tp5,'r*') % fp,tp rate at thresh=0.5
set(gca,'fontsize',12)
title("ROC curve " + test_subj_data.name + " (AUC = "+string(AUC)+")")
xticks(0:0.1:1)
xlabel("False positive rate")
ylabel("True positive rate")

axis square % makes axes equal lengths, so image is undistorted

end



function auc = calc_AUC(false_pos, true_pos)
auc = 0;
for i=2:length(false_pos)
    delta = false_pos(i)-false_pos(i-1);
    auc = auc + 0.5*(true_pos(i-1) + true_pos(i))*delta;
end

end



function [false_pos, true_pos] = calc_fp_tp(thresh, est_probs, test_labels)
% thresh: is the posterior probability threshold (between 0 and 1) above which we
% classify a contact as wm

% est_probs: est_probs{shank_i} is the marginal posterior prob that each contact
% is a wm contact

% test_labels: test_labels{shank_i} gives the true labels for shank_i (1=wm,
% 0=gm)
n_shanks = length(est_probs);
n_fp = 0; n_tp = 0; % number of false pos and true pos
total_p = 0; total_n = 0; % total number of positive and negative
for shank_i=1:n_shanks
    n_fp = n_fp + sum((est_probs{shank_i}(:) > thresh).*(~logical(test_labels{shank_i}(:))));
    n_tp = n_tp + sum((est_probs{shank_i}(:) > thresh).*(logical(test_labels{shank_i}(:))));
    total_p = total_p + sum(test_labels{shank_i});
    total_n = total_n + sum(~logical(test_labels{shank_i}));
end

false_pos = n_fp/total_n;
true_pos = n_tp/total_p;

end
