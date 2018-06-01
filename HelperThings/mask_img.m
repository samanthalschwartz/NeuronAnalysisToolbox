function mask = mask_img(img_in)
thresh = multithresh(single(img_in),2);
mask = img_in>thresh(1);
end