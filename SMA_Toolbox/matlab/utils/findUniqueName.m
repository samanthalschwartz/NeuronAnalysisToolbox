% findUniqueName.m - Safely save a value as a globabl variable
% Mark J. Olah (mjo@cs.unm.edu)
% October 2016
%

function name=findUniqueName(current_names, pattern,idx)
    % Finds a unique name that matches pattern and is distinct from all names in current_names
    % IN:
    %  current_names: cell array of strings of names that are already in use and should not be duplicated
    %  pattern: string containing a single '%i' pattern that represents the pattern of name desired
    %  idx: [optional] starting index [default=1]
    % OUT:
    %  name - type:string - A new unique name matching the pattern with the first unused idx substituted
    if nargin<3
        idx=1;
    end
    idxs = cellmap(@(n) sscanf(n,pattern), current_names);
    idxs = [idxs{:}];
    if ~isempty(idxs)
        idx = max(max(idxs)+1,idx);
    end
    name = sprintf(pattern, idx);
end
