% readmultitif() - A function to read multi-page .tif files
% Mark J. Olah (mjo@cs.unm.edu)
% October 2014
%
function im=readmultitif(fname,waitH)
    % This is like dip_image readtimeseries, but does not read sequentially
    % numbered files, only multi-page tifs.  This gets around bugs in
    % readtimeseries which currently make it impossible to read a mulitpage
    % tif file that ends in a number. Tested to be approx 10-30% slower
    % than readtimeseries as we have to read each frame (page) one at a
    % time.
    %
    % In:
    % fname - type: string - The complete path to the TIFF file including extension
    % waitH - [OPTIONAL] type: handle - A handle to a waitbar ui object to update 
    %
    % Out:
    % im - type:uint16 Size: [Y, X, T] - The image stack, a 3D matlab array
    if nargin==1
        waitH=[];
    end
    
    [~,file_base, ext] = fileparts(fname);
    switch ext
        case {'.tiff','.tif','.TIFF','.TIF'}
            %OK
        otherwise
            error('readmultitif:input','File is not a .tiff "%s"', fname);
    end
    endsInNumber =  '0' <= file_base(end) && file_base(end)<='9';%If it ends in a number readtimeseries won't work
    info=dipio_imagefilegetinfo(fname,'tiff',0);
    nFrames = info.numberOfImages;
    dipioMaxFrames = 4000; %Dipio can't read more than this.

    java_runtime = javaMethod('getRuntime', 'java.lang.Runtime');
    maxJavaMemory = java_runtime.maxMemory() / (1024 * 1024);
    
    tic;
    if nFrames <= dipioMaxFrames
        if endsInNumber
            im=readmultitiff_dipio(fname, waitH);
        else
            im=readmultitiff_readtimeseries(fname, waitH);
        end
    else
        im=readmultitiff_imread(fname, waitH);
%         if exist('bfopen','file')~=2 || maxJavaMemory<=512
%             %do the imread option
%             im=readmultitiff_imread(fname, waitH);
%         else
%             try
%                 im=readmultitiff_bfopen(fname, waitH);
%             catch
%                 im=readmultitiff_imread(fname, waitH);
%             end
%         end
    end
    fprintf('Tiff size X:%i Y:%i #Frames:%i. Load time %.4g s\n',size(im,2), size(im,1), size(im,3), toc());
end 

function im = readmultitiff_readtimeseries(fname, waitH)
    if ishandle(waitH); waitbar(1/100,waitH); end
    im = uint16(readtimeseries(fname));
    if ishandle(waitH); waitbar(100/100,waitH); end
end

function im = readmultitiff_dipio(fname, waitH)
    info=dipio_imagefilegetinfo(fname,'tiff',0);
    N = info.numberOfImages
    %Check for a bad last frame as seems to happen with OHSU tif data
    try
        [~]=dipio_imagereadtiff(fname,N-1);
    catch
        N=N-1;
    end
    im=zeros(info.size(2), info.size(1), N,'uint16');
    for n=1:N %copy each page individually.
        if ishandle(waitH); waitbar(n/N,waitH); end
        im(:,:,n)=uint16(dipio_imagereadtiff(fname,n-1));
    end 
end

function im = readmultitiff_bfopen(fname, waitH)
    data = bfopen(fname);
    data = data{1}; % Only have one "series"    
    N=numel(data);
    im = zeros(size(data,2), size(data,1), N, 'uint16');
    for n=1:N
        if ishandle(waitH); waitbar(n/N,waitH); end
        im(:,:,n) = uint16(data{n});
    end
end

function im = readmultitiff_imread(fname, waitH)
    info = imfinfo(fname);
    N = numel(info);
    im = zeros(info(1).Height, info(1).Width, N,'uint16');
    for n=1:N
        if ishandle(waitH); waitbar(n/N,waitH); end
        im(:,:,n) = uint16(imread(fname, 'Index', n));
    end
end
