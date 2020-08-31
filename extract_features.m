% A function to extract the features needed by the classifier for a given
% patient's data

function out = extract_features(subj_data)

% subj_data: a matlab struct with the following fields:

% name: a string, giving the patient ID or name, i.e. "patient 03"

% data:  < n_timesteps, n_contacts >  array, giving the recorded voltage on
% every contact at each timestep

% ch_type:  < 1, n_contacts > array corresponding to the columns of data,
% with each entry 0, 1, or -1 indicating whether each contact is 
% gray matter (0), white matter (1), or other (-1).
% "other" includes contacts that are outside the brain, malfunctioning, etc

% ch_name:  < 1, n_contacts > cell array, giving the contact name (as a string or char array)
% for each contact (corresponding to the columns of data)

% shank_elecs:  < 1, n_shanks > cell array, where shank_elecs{i} is an array
% of contact indices for shank i. These are the columns of data that correspond
% to the contacts on a particular shank. THESE MUST BE ORDERED FROM SHALLOW
% TO DEEP; i.e. shank_elecs{1}(3) is shallower in the brain (closer to
% skull) than shank_elecs{1}(4).

% Fs:  the sampling rate used to record this patient's data (in Hz)


% make initial ch_name_list, ch_type_list
ch_name = subj_data.ch_name;
ch_type = subj_data.ch_type;
shank_elecs = subj_data.shank_elecs;

ch_name_list = {};
ch_type_list = {}; 
n_shanks = length(shank_elecs);
for shank_i=1:n_shanks
    ch_name_list{shank_i} = {ch_name{shank_elecs{shank_i}}}';
    ch_type_list{shank_i} = ch_type(shank_elecs{shank_i})';
end


% make bipolar data and calc feature (ch_feat_list)
% In both EFRI and LA data, the highest numbered channel (i.e. A10, B10, etc) on a
% shank is the outermost channel (closest to skull), and therefore is often
% 'out'. In EFRI, this is the channel at the beginning of each shank's list (index
% 1) (in LA it's the last channel, index end), so here we start at index 1 and subtract off the
% next channel, so channels further in aren't affected by channels further
% out. The first channel is only removed if it's actually "out".
lfpdata = subj_data.data;
Fs = subj_data.Fs;

ch_spec_list_bp = {}; % power spectrum for each channel (bipolar ref)
all_specs_bp = []; % for conveniently computing mean power spectrum (bipolar ref)
n_shanks = length(shank_elecs);
for shank_i=1:n_shanks
    n_ch = length(shank_elecs{shank_i});
    for ch_i=1:n_ch
        if ch_type_list{shank_i}(ch_i) < 0
            ch_spec_list_bp{shank_i}{ch_i} = NaN;
        elseif ch_i == n_ch
            data_bp = lfpdata(:,shank_elecs{shank_i}(ch_i)) - ...
                lfpdata(:,shank_elecs{shank_i}(ch_i-1));
            ch_spec_list_bp{shank_i}{ch_i} = calc_spectrum(data_bp, Fs);
            all_specs_bp(end+1,:) = ch_spec_list_bp{shank_i}{ch_i};
            
        else
            data_bp = lfpdata(:,shank_elecs{shank_i}(ch_i)) - ...
                lfpdata(:,shank_elecs{shank_i}(ch_i+1));
            ch_spec_list_bp{shank_i}{ch_i} = calc_spectrum(data_bp, Fs);
            all_specs_bp(end+1,:) = ch_spec_list_bp{shank_i}{ch_i};
            
        end
    end
end

mean_spectrum_bp = mean(all_specs_bp,1);
% remove 4 hz band around 60 and 120 hz
mask_bp = logical(ones(1,length(mean_spectrum_bp))); mask_bp(58:62)=false; mask_bp(118:122)=false;
ch_feat_list = {};
n_shanks = length(shank_elecs);
for shank_i=1:n_shanks
    outermost_ch = find(ch_type_list{shank_i} > -1, 1); % index of outermost contact
    ch_feat_list{shank_i} = [];
    n_ch = length(shank_elecs{shank_i});
    for ch_i=1:n_ch
        if ch_type_list{shank_i}(ch_i) < 0
            ch_feat_list{shank_i}(ch_i,:) = [NaN, NaN];
        else
            ch_feat_list{shank_i}(ch_i,:) = [mean(ch_spec_list_bp{shank_i}{ch_i}(mask_bp) - mean_spectrum_bp(mask_bp)),...
                ch_i - outermost_ch];
        end
    end
end


% remove invalid channels: channels with type =-1
n_shanks = length(shank_elecs);
for shank_i=1:n_shanks
    n_ch = length(shank_elecs{shank_i});
    for ch_i=n_ch:-1:1 % need to go in reverse order bc we're changing the number of channels as we go
        if ch_type_list{shank_i}(ch_i) < 0
            ch_name_list{shank_i}(ch_i) = [];
            ch_type_list{shank_i}(ch_i) = [];
            ch_feat_list{shank_i}(ch_i,:) = [];
        end
    end
end


% remove invalid shanks: shanks with 0 channels
n_shanks = length(shank_elecs);
for shank_i=n_shanks:-1:1
    if length(ch_name_list{shank_i}) == 0
        ch_name_list(shank_i) = []; % use round brackets when deleting elements of cell array
        ch_type_list(shank_i) = [];
        ch_feat_list(shank_i) = [];
    end
end


out.name = subj_data.name;
out.ch_name_list = ch_name_list;
out.ch_type_list = ch_type_list;
out.ch_feat_list = ch_feat_list;

end




function spectrum = calc_spectrum(data_bp, Fs)
% data_bp is assumed to be a 1D array of lfpdata for 1 electrode
n_samples = 10;
sample_length = 10; % in sec

samples = {};
last_indx = round(length(data_bp)-sample_length*Fs);
start_indices = round(linspace(1,last_indx,n_samples));
stop_indices = start_indices + round(sample_length*Fs) - 1;
for i=1:n_samples
    samples{i} = data_bp(start_indices(i):stop_indices(i));
end

freqs = 1:1:150;

all_pxx = nan(n_samples,length(freqs));
for i=1:n_samples
    [all_pxx(i,:), f] = periodogram(samples{i},[],freqs,Fs);
end

spectrum = mean(log10(all_pxx),1);
end







