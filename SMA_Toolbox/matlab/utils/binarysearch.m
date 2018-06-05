function i=binarysearch(svec, x)
    % [IN]
    % svec - sorted vector
    % x - value to look for
    % [Out]
    % i - index of insert 1 more than the last index that is strictly smaller than x
    
    a=1;
    b=numel(svec);
    if x <=svec(a)
        i=1;
    elseif x > svec(b)
        i = b+1;
    else
        while a+1<b        
            c = floor((a+b)/2);
            if x <= svec(c)
                b = c;
            else
                a = c;
            end
        end
        i = a+1;
    end
end