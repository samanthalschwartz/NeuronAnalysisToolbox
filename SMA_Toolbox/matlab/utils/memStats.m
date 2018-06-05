

function [mem] = memStats()
    if ~isunix
    mem.error = true;
    return
    end

    [sts,msg] = unix('free -m | grep Mem:');

    if sts %error
        mem.error = true;
        return
    else %no error
    mems = cell2mat(textscan(msg,'%*s %u %u %u %*u %*u %*u','delimiter',' ','collectoutput',true,'multipleDelimsAsOne',true));
    mem.freeMB = mems(3);
    mem.usedMB = mems(2);
    mem.totalMB = mems(1);
    mem.error = false;
end
