# WM-classifier
Classification of stereo-EEG contacts in white matter versus gray matter using recorded activity

A description of the method can be found in the accompanying paper.

To use, each patient's data must be saved as a separate matlab struct with the following fields:

name: a string, giving the patient ID or name, i.e. "patient 03"

data:  < n_timesteps, n_contacts >  array, giving the recorded voltage on
every contact at each timestep

ch_type:  < 1, n_contacts > array corresponding to the columns of data,
with each entry 0, 1, or -1 indicating whether each contact is 
gray matter (0), white matter (1), or other (-1).
"other" includes contacts that are outside the brain, malfunctioning, etc

ch_name:  < 1, n_contacts > cell array, giving the contact name (as a string or char array)
for each contact (corresponding to the columns of data)

shank_elecs:  < 1, n_shanks > cell array, where shank_elecs{i} is an array
of contact indices for shank i. These are the columns of data that correspond
to the contacts on a particular shank. THESE MUST BE ORDERED FROM SHALLOW
TO DEEP; i.e. shank_elecs{i}(3) is shallower in the brain (closer to
skull) than shank_elecs{i}(4).

Fs:  the sampling rate used to record this patient's data (in Hz)


Features can then be extracted and saved by running extract_features_script.m
The classifier is then run on a given patient using main_script.m
