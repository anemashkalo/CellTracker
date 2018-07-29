function [newstr,raw_im,nuc_dat1,nuc_dat2]=get_celltypes_centroids(all_h5_type2,mask_out2,binary_mask1,seprate_mask2_clean)

raw_im = struct;
nuc_dat1 = struct;
nuc_dat2 = struct;
stats1 = [];
stats2 = [];

for ii =1:size(mask_out2,2)
newstr = [all_h5_type2(ii).file(1:end-17) '.tif' ]; % dapi image
raw_im(ii).dat = imread(newstr);
% get centroids of two cell types
stats1 = regionprops(binary_mask1(ii).dat,'Centroid','Area');
stats2 = regionprops(seprate_mask2_clean(ii).dat,'Centroid','Area');

nuc_dat1(ii).xy = cat(1,stats1.Centroid);
nuc_dat1(ii).area = cat(1,stats1.Area);

nuc_dat2(ii).xy = cat(1,stats2.Centroid);
nuc_dat2(ii).area = cat(1,stats2.Area);



end

end