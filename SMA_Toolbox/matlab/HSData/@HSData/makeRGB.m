function im=makeRGB(obj, cm)
    %
    % Tested at 20x speed of HSIData.makeRGB(), and with significantly cleaner
    % background.  Any RGB value less than 0 is trucated to 0.
    %
    % [in] cm - a Nx3 array of colors mapping the Nth wavelength to the Nth color.
    % [out] im - [y,x] or [y,x,t] RGB dip_image or dip_image sequence.
    %
    if nargin<2
        cm = obj.colorMap;
    end
    %Make correct sized channel arrays
    spatialSize = [obj.sizeX,obj.sizeY,obj.nFrames];
    pixels = reshape(obj.getFrames(), obj.sizeL, prod(spatialSize)); % treat as 2D matrix of row:L  col:Pixel

    R = zeros(spatialSize);
    G = zeros(spatialSize);
    B = zeros(spatialSize);

    % Eliminate negative pixels
    R(:) = max(0, cm(:,1)' * pixels);
    G(:) = max(0, cm(:,2)' * pixels);
    B(:) = max(0, cm(:,3)' * pixels);

    %Normalize
    R = R/max(R(:));
    G = G/max(G(:));
    B = B/max(B(:));
    
    im = joinchannels('RGB', R, G, B);
end
