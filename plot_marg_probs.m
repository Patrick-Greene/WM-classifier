% plot marginal probabilities
function out = plot_marg_probs(est_probs, test_subj_data)

colormap(brewermap(100,'*RdBu'));
clims = [0,1]; % max min vals for colormap
n_shanks = length(est_probs);
j=1;
for shank_i=1:n_shanks
    if length(est_probs{shank_i}) > 0 % skip shanks w/ no elecs
        subplot(n_shanks,2,j)
        imagesc(est_probs{shank_i}, clims)
        shank_elecs = test_subj_data{1}.ch_name_list;
        xticks(1:length(shank_elecs{shank_i}))
        xticklabels(shank_elecs{shank_i})
        yticklabels([])
        if j==1
            title("Probability of white matter",'FontSize',14,'FontWeight','Normal')
        end
        j = j+1;
        subplot(n_shanks,2,j)
        % lighten the colors slightly on the true labels so easier to see
        imagesc((0.9*double(test_subj_data{1}.ch_type_list{shank_i}'))+0.05, clims) 
        xticks(1:length(shank_elecs{shank_i}))
        xticklabels(shank_elecs{shank_i})
        yticklabels([])
        if j==2
            title("True label",'FontSize',14,'FontWeight','Normal')
        end
        j=j+1;
    end
end
% sgtitle("Posterior probability subject " + all_subj_data{subj_i}.subj_data.name)
sgtitle(test_subj_data{1}.name,'FontSize',14,'FontWeight','Normal')

end