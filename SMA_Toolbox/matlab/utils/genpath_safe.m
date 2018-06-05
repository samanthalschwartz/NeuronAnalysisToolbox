% saveGlobalVar.m - Generate recursive path strings that exclude patterns
% Mark J. Olah (mjo@cs.unm.edu)
% October 2014
%
function p=genpath_safe(rootdir, exclude_re)
    % Like genpath but allows the exclusion of directories matching any of
    % the reguluar expressions in exclude_re.
    %
    % This should work on windows and linux, we deal with the differnces in
    % matlab path string format as well as with the directory seperator
    % format.
    %
    % IN:
    %  rootdir - type:string - The directory to start the search at
    %  exclude_re - [OPTIONAL] cell array of strings - The regular
    %               expressions to exclude
    % OUT:
    % p - type:string - The path string suitable for addpath.
    
    if ispc()
        spch=';';
    else
        spch=':';
    end
    if nargin==1
        exclude_re={'[\\/]\.svn'};
    end
    p=genpath(rootdir);
    if ~isempty(p)
        ps=strsplit(p,spch);
        for n=1:length(exclude_re)
            re=exclude_re{n};
            m=regexp(ps, re, 'match');
            ps=ps(cellfun(@isempty,m));
            if isempty(ps);
                p=[];
                return
            end
        end
        p=strjoin(ps,spch);
    end
end