% Mark J. Olah (mjo@cs.unm.edu)
% 02/2015


function rpath=relativepath(source_dir, target)
    % [IN]
    % source_dir - full path to source directory
    % target - full path to target (must trail with a slash or a file name
    % [OUT]
    % rpath - the relative path using the ".." notation including the listed file from the target
    [target_path,file,ext]=fileparts(target);

    S=strsplit(source_dir,filesep());
    S(cellfun(@isempty,S))=[]; % remove empty cells
    T=strsplit(target_path,filesep());
    T(cellfun(@isempty,T))=[]; % remove empty cells
    N=min(length(S),length(T));
    n=find(~cellfun(@strcmp,S(1:N),T(1:N)),1,'first');
    if isempty(n)
        n=N+1;
    end
    R=cell(1,length(S)-n+1);
    R(:)={'..'};
    R=[R T(n:end)];
    rpath=collapsepath(fullfile(R{:},[file ext]));
end
