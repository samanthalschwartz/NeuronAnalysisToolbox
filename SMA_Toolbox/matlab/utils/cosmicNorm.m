function normIm = cosmicNorm(im, k)
    % Normalization for data subject to cosmic rays.  This normalizes the image to have a maxima of 1 at
    % the k-th largest value.  All values larger than the maxima are also set to 1.
    if nargin==1
        k=16;
    end
    normIm = zeros(size(im));
    normIm(:) = max(0, im);
    normIm = normIm - min(normIm(:));
    sk=sort(normIm(:), 'descend');
    mx=sk(k);
    normIm(:) = min(normIm(:),mx)./mx;
end
