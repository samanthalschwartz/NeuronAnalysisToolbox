

path = '/home/mjo/LidkeLab/Data/Emanuel/150619_A431_EGF-QD585_EGF-QD655/';
spdatas = {'cell1_160pM-EGF-QD-each_set3_500fps#0001-2015-6-19-15-53-58.spdata'};

full_spdatas = cellmap(@(f) fullfile(path,f), spdatas); %These are the full paths to spdatas





data_dir=fullfile(getenv('LIDKE_LAB'), 'Data', 'HSM');
spthsi_file=fullfile(data_dir,'test.spthsi'); %This should ba a SPTHSI file

warning('OFF', 'SPTloadFile:dataFileNotLinked'); %Don't need linked data
spthsi=SPTHSI(spthsi_file); %Make an SPTHSI object from file
tgrp=HSTrackGroup(spthsi); %Make an HSTrackGroup from SPTHSI 

%Set up general parameters
default_certainty=0.95 %default certatinty to use when testing distances with known stddev

%Set up track filtering critera
track_min_length=10; %min track localization count
track_min_duration=1.0; %min track duration (s)
track_min_duty_ratio=0.25; %min track duty ratio (ratio of localized to total frames)
track_max_spectral_dev=15; %max specral deviation (nm)

%Set up pair filtering criteria
pair_min_length=10; %min pair simultaneous (same frame) localization count
pair_min_duration=1.0; %min pair simultaneous duration (s) [i.e. temporal overlap]
pair_min_spectra_dist

%Filter tracks
fprintf('Using SPTHSI file: %s\n', spthsi_file);
fprintf('HSTrackGroup: #Tracks:%i #Localizations:%i\n',tgrp.nT, tgrp.nL);
% fprintf('HSTrackGroup: start_time:%g end_time:%g duration:%g\n',tgrp.nT, tgrp.nL);

nrejected=tgrp.filter_track_length(track_min_length);
fprintf('Rejected %i tracks with length<%i\n',nrejected,track_min_length);

nrejected=tgrp.filter_track_duration(track_min_duration);
fprintf('Rejected %i tracks with duration<%.3fs\n',nrejected,track_min_duration);

nrejected=tgrp.filter_track_duty(track_min_duty_ratio);
fprintf('Rejected %i tracks with duty ratio<%.2f\n',nrejected,track_min_duty_ratio);

nrejected=tgrp.filter_track_spectral_deviation(track_max_spectral_dev, certatinty);
fprintf('Rejected %i tracks with spectral deviation>%.1fnm (@%.1f%% conf.).\n',nrejected,track_max_spectral_dev, certainty*100.);
fprintf('Filtered HSTrackGroup: #Tracks:%i #Localizations:%i\n\n',tgrp.nT, tgrp.nL);

%Make all pairs that overlap in time
pairs=tgrp.make_temporal_pairs(pair_min_duration);
npairs=size(pairs,1);
fprintf('Found %i pairs with min duration<%.3fs',npairs,min_pair_duration);

pairs=tgrp.filter_pair_length(pairs,pair_min_length);
npairs=size(pairs,1);
fprintf('Rejected %i pairs with min length<%i',size(pairs,1),min_pair_duration);

pairs=tgrp.filter_pair_spectra_dist(pairs,pair_min_spectra_dist, certainty);

pairs=tgrp.filter_pair_spectra_dist(pairs,pair_min_spectra_dist, certainty);

tgrp.print();










