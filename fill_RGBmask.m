function out_mask = fill_RGBmask(mask, mask_2D,RGB_val)


for i = 1:3
    temp(:,:,i) = ones(size(mask,1),size(mask,2))*RGB_val(i);
end

mask(repmat(mask_2D,[1,1,3])) = temp(repmat(mask_2D,[1,1,3]));
out_mask = mask;
