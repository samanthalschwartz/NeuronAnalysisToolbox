%--- example file to add the correct paths

dippath = '/Applications/dip/common/dipimage';
addpath(dippath);
dip_initialise
dipsetpref('DefaultGlobalStretch',true);

path2ForChris = 'put path here';
addpath(path2ForChris);

path2loadtiff =  'put path here';
addpath(path2loadtiff);

path2startupSMA = 'put path here';
addpath(genpath(path2startupSMA));