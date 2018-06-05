% saveGlobalVar.m - Safely save a value as a globabl variable
% Mark J. Olah (mjo@cs.unm.edu)
% October 2014
%

function name=saveGlobalVar(name,val)
    %Saves a value as global variable name, but adds an index number to
    %prevent overwriting current global variables
    % IN:
    %  name - type:string - The globabl variable name to use
    %  val - any value
    % OUT:
    %  name - type:string - The actual name used with index added
    cmd = sprintf('who(''%s*'')',name);
    globals = evalin('base',cmd);
    if ~isempty(globals)
        re = sprintf('%s(?<idx>\\d+)',name);
        [idxs,~] = regexp(globals,re,'tokens','match');
        idxs = idxs(~cellfun(@isempty,idxs));
        if isempty(idxs)
            idx = 1;
        else
            idx = 1+max(cellfun(@(a) str2double(a{1}),idxs));
        end
        name = sprintf('%s%i',name,idx);
    end
    assignin('base',name,val);
end
