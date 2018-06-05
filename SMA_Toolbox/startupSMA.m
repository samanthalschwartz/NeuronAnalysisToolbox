%
% Mark J. Olah (mjo@cs.unm.edu)
% 09/2014
%

function startupSMA( debug, sma_path)
    %
    % Sets up the paths to matlab and mex code for the SMA_Toolbox for windows and linux
    %
    % Inputs:
    %   debug - boolean, determines if we should use debugging version of mex libararies or not (Default: false)
    %   sma_path - The full path to the SMA_Toolbox folder.  Defaults to the location of this setupSMA file.
    %

    if nargin==0
        debug=false;
    end
    if nargin<=1
        sma_path=strsplit(which('startupSMA'),'%');% Remove comments from which()
        [sma_path,~,~]=fileparts(sma_path{1});
    end

    if ispc()
        if debug
            mex_sub_dir='mex.w64.debug';
        else
            mex_sub_dir='mex.w64';
        end
    elseif isunix()
        if debug
            mex_sub_dir='mex.glnxa64.debug';
        else
            mex_sub_dir='mex.glnxa64';
        end
    else
        error('Platform not supported.');
    end
    sma_mex_path=fullfile(sma_path,'mex',mex_sub_dir);
    sma_matlab_path=fullfile(sma_path,'matlab');
    addpath(sma_mex_path);
    addpath(fullfile(sma_matlab_path, 'utils')); %this is where genpath_safe is
    addpath(genpath_safe(sma_matlab_path));
end
