function sv_arr2 = remakeRegistration_forChannelImagedEveryNFrames(original_sv_arr,numframes2skipAtBeginning,increment)
sv_arr2 = zeros(size(original_sv_arr));
rem = mod(size(original_sv_arr,2)+1-numframes2skipAtBeginning,increment);
% reshape sv_arr for ch2:
for tt = (numframes2skipAtBeginning):increment:(size(original_sv_arr,2)-rem)
    sv_arr2(:,tt:(tt+(increment-1))) = repmat(original_sv_arr(:,tt),[1 increment]);    
end
sv_arr2(:,(size(original_sv_arr,2)-rem+1):end) = repmat(original_sv_arr(:,size(original_sv_arr,2)-rem+1),[1 rem]);