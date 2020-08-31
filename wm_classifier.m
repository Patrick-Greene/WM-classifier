function est_probs = wm_classifier(train_features_wm, train_features_gm,...
    train_labels, test_features, max_ch, alpha_cd, n_samples)
% This function integrates over the posterior distribution of the parameters
% to get the posterior predictive

% train_features_gm/wm: training features for every contact on every training patient, concatenated
% together. Each row is a contact, each column a feature.

% test_features: test_features{shank_i} are the features for shank_i; each
% row is a contact, each column a feature

% prior_pdf: a function that, given a set of labels for all contacts on a
% shank, returns the value of the prior pdf

% alpha_cd : kernel bandwidth for ch dist feature (user selected)
% n_samples : number of samples to use when integrating over param
% posterior

% Returns:
% est_probs: est_probs{shank_i} is the marginal posterior predictive prob that each contact
% is a wm contact


% ch_dist_indx is the index of the feature component which is the integer distance
% from top of shank. I.e. if this is the 2nd feature component, then
% ch_dist_indx=2). Set disc_indx=0 if not using this as a feature or don't want to treat
% it as discrete.
ch_dist_indx = 2;

% Calc peak and covariance of parameter posterior dist p(alpha,beta|D)
alpha_wm0 = 0.2; % initial val for optimization
alpha_gm0 = 0.2;
beta0 = 3;
n_feat = length(alpha_wm0);
[mu, sigma] = calc_max_params(alpha_wm0, alpha_gm0, beta0,...
    train_features_wm,train_features_gm, train_labels, max_ch, ch_dist_indx, alpha_cd);

% get samples
samples = gauss_approx_posterior_params_sampler(n_samples, mu, sigma);
% each row is a sample of the form [alpha_wm alpha_gm beta]
% Get normalizing const for each beta sample
sample_norm_const = {};
for sample_i=1:n_samples
    sample_norm_const{sample_i} = calc_normalizing_constants(max_ch, samples(sample_i,end));
end

est_probs = {}; % est_probs{shank_i} is a vector of marginal probabilities that each electrode
% is white matter
n_shanks = length(test_features);
for shank_i=1:n_shanks
    marg_prob_samples = [];
    for sample_i=1:n_samples
        marg_prob_samples(end+1,:) = marg_shank_posterior(train_features_wm, train_features_gm,...
            test_features{shank_i}, ch_dist_indx, @prior_pdf5, alpha_cd,...
            samples(sample_i,1:n_feat), samples(sample_i,n_feat+1:2*n_feat), samples(sample_i,end),...
            sample_norm_const{sample_i});
    end
    est_probs{shank_i} = mean(marg_prob_samples,1);
end

end



function [mu, sigma] = calc_max_params(alpha_wm0, alpha_gm0, beta0,...
    features_wm,features_gm, labels, max_ch, ch_dist_indx, alpha_cd)
% calc peak and covariance of p(alpha,beta|D).
% mu: [mu_alpha_wm mu_alpha_gm mu_beta]

% first optimize over alpha_wm and alpha_gm separately, since independent
n_feat = length(alpha_wm0);
options = optimoptions('fmincon','StepTolerance',0.001,'Display','off');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) -log_un_posterior_alpha_generic(...
    x,features_wm, ch_dist_indx, alpha_cd),...
    alpha_wm0, -eye(n_feat), zeros(n_feat,1), [],[],[],[],[], options);

mu = x;
sigma = zeros(2*n_feat+1);
sigma(1:n_feat,1:n_feat) = inv(hessian); % If 1D, hessian is just a number

options = optimoptions('fmincon','StepTolerance',0.001,'Display','off');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) -log_un_posterior_alpha_generic(...
    x,features_gm, ch_dist_indx, alpha_cd),...
    alpha_gm0, -eye(n_feat), zeros(n_feat,1), [],[],[],[],[], options);

mu = [mu x];
sigma(n_feat+1:2*n_feat,n_feat+1:2*n_feat) = inv(hessian);

% optimize over beta
options = optimoptions('fmincon','StepTolerance',0.001,'Display','off');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) -log_un_posterior_beta(...
    x, labels, max_ch), [beta0], -eye(1), zeros(1), [],[],[],[],[], options);

mu = [mu x];
sigma(end,end) = inv(hessian);
    
end


function est_probs = marg_shank_posterior(train_features_wm, train_features_gm,...
    test_features, ch_dist_indx, prior_pdf, alpha_cd, alpha_wm, alpha_gm, beta, norm_const)

% train_features_gm/wm: features for every contact on every patient, concatenated
% together. Each row is a contact, each column a feature. Since just one
% feature right now, it's just a long column vector.

% test_features: 1 x n_feats the features for 1 shank

% ch_dist_indx: the index of the feature component which is the integer distance
% from top of shank. I.e. if this is the 3rd feature component, then
% disc_indx=3). Set disc_indx=0 if not using this as a feature or don't want to treat
% it as discrete.

% prior_pdf: a function that given a set of labels for all contacts on a
% shank, returns the value of the prior pdf. Must also supply beta and
% corresponding normalizing constants.

% Returns:
% est_probs: est_probs is the marginal posterior prob that each contact
% is a wm contact for the given shank


% Classify all the electrodes on each shank of the excluded subject
% est_probs is a vector of marginal probabilities that each electrode
% is white matter
% gm_lik is likelihood for all elecs on shank i, under the gm
% distribution. Similarly for wm_lik but under wm distribution

gm_lik = calc_likelihood2(train_features_gm, test_features, ch_dist_indx, alpha_gm, alpha_cd);
wm_lik = calc_likelihood2(train_features_wm, test_features, ch_dist_indx, alpha_wm, alpha_cd);

% generate all possible wm/gm labelings of the electrodes on shank_i,
% and the associated likelihood*prior (un_post) for each
[labelings, un_post] = gen_labelings2(gm_lik, wm_lik,@(x) prior_pdf(x, beta, norm_const));
normalizing_const = sum(un_post); % bayes normalizing constant

% get the marginal probability estimate
n_elecs = length(test_features);
marg_probs = nan(1,n_elecs);
for elec_i=1:n_elecs
    % Note labelings(:,elec_i) gives all the rows where elec_i == 1
    marg_probs(elec_i) = sum(un_post(logical(labelings(:,elec_i))));
end
est_probs = marg_probs/normalizing_const;
      
end


function samples = gauss_approx_posterior_params_sampler(n_samples, mu, sigma)
% mu: [mu_alpha_wm mu_alpha_gm mu_beta] and covariance for truncated gaussian approx to posterior param
% distribution. Can be used directly or as part of importance sampler
% returns: n_samples x 2*n_continuous_feats + 1 array. I.e. if we have bp
% and cr as continuous features, then we're drawing samples from a 5D
% distribution, since we have to get means for both the wm and gm kernel
% widths

samples = nan(n_samples,length(mu));
for i=1:n_samples
    is_pos = false; % check if all components are positive
    while is_pos == false
        x = mvnrnd(mu, sigma);
        if all(x > 0)
            samples(i,:) = x;
            is_pos = true;
        end
    end
end
end


function out = log_un_posterior_alpha_generic(alpha,features, ch_dist_indx, alpha_cd)
lambda = 0.1; % This is also defined in log_prior_alpha (make sure they match)

loglik_alpha = calc_alpha_loglik_generic(features, ch_dist_indx, alpha, alpha_cd);
out = loglik_alpha + log_prior_alpha_generic(alpha, lambda);
end

function out = log_un_posterior_alpha(alpha_wm, alpha_gm,features_wm,...
    features_gm, ch_dist_indx, alpha_cd)
loglik_alpha = calc_alpha_loglik(features_wm, features_gm, ch_dist_indx, alpha_wm, alpha_gm, alpha_cd);
out = loglik_alpha + log_prior_alpha(alpha_wm, alpha_gm);
end

function out = log_un_posterior_beta(beta, labels, max_ch)
loglik_beta = calc_beta_loglik(beta, labels, max_ch);
out = loglik_beta + log_prior_beta(beta);
end


function out = log_un_posterior_params(alpha_wm, alpha_gm, beta, features_wm,...
    features_gm, labels, max_ch, ch_dist_indx, alpha_cd)
% unnormalized log posterior distribution of alpha, beta params
% features_wm/gm: N x n_feats array concatenating the
% feature values for all patients and all contacts that are either gm or wm
% labels{i} is binary vector of labels for some shank i of some subject in training set
loglik_alpha = calc_alpha_loglik(features_wm, features_gm, ch_dist_indx, alpha_wm, alpha_gm, alpha_cd);
loglik_beta = calc_beta_loglik(beta, labels, max_ch);

out = loglik_alpha + loglik_beta + log_prior_alpha(alpha_wm, alpha_gm) + log_prior_beta(beta);
end


function out = log_prior_alpha_generic(alpha, lambda)
out = 0;
for i=1:length(alpha)
    out = out + -lambda*alpha(i);
end
end


function out = log_prior_alpha(alpha_wm, alpha_gm)
% unnormalized log prior dist on alpha, the kernel widths for the continuous features
lambda = 0.1;
out = log_prior_alpha_generic(alpha_wm, lambda) + log_prior_alpha_generic(alpha_gm, lambda);
end

function out = log_prior_beta(beta)
% unnormalized log prior dist on beta, the smoothing parameter for the prior
lambda = 0.1;
out = -lambda*beta;
end


function loglik = calc_alpha_loglik_generic(features, ch_dist_indx, alpha, alpha_cd)
% features: either features_wm or features_gm
% alpha: either alpha_wm or gm

n_data = size(features,1);

loglik = 0;
for i=1:n_data
    mask = false(n_data,1); mask(i)=true;
    loglik = loglik + log(calc_likelihood2(features(~mask,:), features(mask,:),...
        ch_dist_indx, alpha, alpha_cd));
end

end

function loglik = calc_alpha_loglik(features_wm, features_gm, ch_dist_indx, alpha_wm, alpha_gm, alpha_cd)
% log p(x_train | z_train, alpha) where alpha is vector of kernel widths
% alpha_gm(i), alpha_wm(i) are kernel widths for feature i of gm/wm dist

loglik = calc_alpha_loglik_generic(features_wm, ch_dist_indx, alpha_wm, alpha_cd) + ...
    calc_alpha_loglik_generic(features_gm, ch_dist_indx, alpha_gm, alpha_cd);
end


function lik = calc_likelihood2(train_features, test_features, ch_dist_indx, alpha, alpha_cd)
% Same as calc_likelihood, but allows kernel width alpha to be specified
% rather than using silverman's rule of thumb.
% train_features: N x n_features array, giving features for all wm or gm training contacts
% test_features: n_ch x n_features feature array, giving features for each
% contact on a particular shank
% ch_dist_indx: the index in the feature vector of the component representing distance from top of shank
% that we should integrate over. I.e. ch_dist_indx=3 means the third
% feature component is the channel distance. If ch_dist_indx == 0, it's either not a 
% feature or we're not integrating over it.
% alpha: scalar or row vector of kernel widths for the continuous
% feature(s). If train_features is train_features_wm, then alpha is
% alpha_wm, and same with gm.
% alpha_cd: ch distance kernel width (user specified)

% lik: n_ch x 1 vector of likelihoods

% Since distance from top of shank is really an integer-valued variable, need to
% integrate from distance-0.5 to distance+0.5

[n_data, n_feat] = size(train_features);
mask = true(1,n_feat); mask(ch_dist_indx)=false;
bandwidth = zeros(1,n_feat); bandwidth(mask) = alpha;
bandwidth(ch_dist_indx) = alpha_cd; % ch distance kernel width

if ch_dist_indx == 0
    % just return values of density without integrating
    lik = mvksdensity(train_features, test_features,'Bandwidth',bandwidth);
else
    n_ch = size(test_features,1);
    lik = zeros(n_ch,1);
    % integrate
    delta_x = 0.25;
    for ch_i=1:n_ch
        ch_dist = test_features(ch_i,ch_dist_indx); % channel distance from top of shank
        if ch_dist == 0
            x = [-3:delta_x:0.5];
        else
            x = [ch_dist-0.5:delta_x:ch_dist+0.5];
        end
        eval_pts = ones(length(x),n_feat).*test_features(ch_i,:);
        eval_pts(:,ch_dist_indx) = x';
        f = mvksdensity(train_features, eval_pts,'Bandwidth',bandwidth);
        lik(ch_i) = sum(0.5*(f(2:end)+f(1:end-1))*delta_x);
    end
end

end 
    

function beta = find_best_beta(all_subj_data, excl_subj_i, max_ch)
% find best beta from training set by maximizing p(beta|z)
% max_ch: max number of channels any shank

labels = {};
n_subj = length(all_subj_data);
for subj_i=1:n_subj
    if subj_i == excl_subj_i
        continue
    end
    n_shanks = length(all_subj_data{subj_i}.subj_data.ch_type_list);
    for shank_i=1:n_shanks
        labels{end+1} = all_subj_data{subj_i}.subj_data.ch_type_list{shank_i};
    end
end

beta = fminbnd(@(x) -calc_beta_loglik(x,labels,max_ch), 0, 12);
end


function loglik = calc_beta_loglik(beta, labels, max_ch)
% calc likelihood p(z|beta) for training set data z
% labels{i} is binary vector of labels for some shank i of some subject in training set
norm_const = calc_normalizing_constants(max_ch, beta);

loglik = 0;
for i=1:length(labels)
    loglik = loglik + log(prior_pdf5(labels{i}, beta, norm_const));
end 
end



function pdf_val = prior_pdf5(labels, beta, norm_const)
% labels: binary labels for a given electrode shank
% norm_const[n] is the normalizing constant for the prior when applied to
% electrodes with n contacts

n_ch = length(labels);
z = 2*labels - 1; % makes z_i = -1,1 instead of 0,1 
% pdf_val = (1/norm_const(n_ch))*exp(beta*sum(z(1:end-1).*z(2:end)));
pdf_val = (1/norm_const(n_ch))*exp((1/n_ch)*beta*sum(z(1:end-1).*z(2:end))); % normalized by 1/n_ch

end


function norm_const = calc_normalizing_constants(max_ch, beta)
% calculate the normalizing constants for electrodes with up to max_ch
% contacts. norm_const[n] is for an n contact shank.
% THESE ARE A FUNCTION OF BETA! So if change beta, need to recalculate!

norm_const = nan(max_ch,1);
for n_ch=1:max_ch
    labelings = dec2bin(0:(2^n_ch)-1)-'0';
    c=0;
    for i=1:2^n_ch
        z = 2*labelings(i,:) - 1; % makes z_i = -1,1 instead of 0,1 
%         c = c + exp(beta*sum(z(1:end-1).*z(2:end)));
        c = c + exp((1/n_ch)*beta*sum(z(1:end-1).*z(2:end))); % normalized by 1/n_ch
    end
    norm_const(n_ch) = c;
end
end


function [labelings, u_post] = gen_labelings2(shank_gm_lik, shank_wm_lik, prior_pdf)
% Same as original, but does calculations on log posterior since the unormalized posterior
% values are extremely small (around 1e-24 - 1e-15), and subtracts off the
% largest (least negative) log posterior value
% shank_gm/wm_lik = gm/wm_lik{shank_i} for whatever shank_i we're doing
% labelings: binary array of all possible labelings of the electrodes (each
% row is a labeling), with wm=1, gm=0
% u_post(i): likelihood*prior (unnormalized posterior) corresponding to row
% i of labelings

n_elecs = length(shank_gm_lik);
labelings = dec2bin(0:(2^n_elecs)-1)-'0';

u_logpost = nan(2^n_elecs,1);
for i=1:2^n_elecs
    loglik = sum(log((shank_wm_lik(logical(labelings(i,:)))))) + ...
        sum(log(shank_gm_lik(~logical(labelings(i,:)))));
    
    logprior = log(prior_pdf(labelings(i,:)));
    
    u_logpost(i) = loglik + logprior;
end
u_logpost = u_logpost - max(u_logpost);
u_post = exp(u_logpost);
end


