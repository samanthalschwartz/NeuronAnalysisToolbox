function times=timeit(count, funs, varargin)
    burnin=1;
    ns=cellfun(@numel,varargin);
    n=ns(1);
    if ~all(ns==n)
        error('timeit:InputSizes','Inputs different lengths');
    end

    if ~iscell(funs)
        funs=cellconst(funs,n); %make n copys of the functions
    elseif n==1 && numel(funs)>1
        n=numel(funs);
        varargin=cellmap(@(inp) cellconst(inp,n), varargin); %make n copys of the args
    end
    inputs=cellmap(@(varargin) varargin(:), varargin{:}); %Think on this one a little.
    %burnin
    for i=1:burnin
        funs{1}(inputs{1}{:});
    end
    %timing
    T=zeros(count,n);
    for c=1:count
        for i=1:n
            tic;
            funs{i}(inputs{i}{:});
            T(c,i)=toc;
        end
    end
    if count>1
        out.mean=mean(T)';
        out.var=var(T)';
        out.std=std(T)';
    else
        out.times=T';
    end
    times=struct2table(out);
end
