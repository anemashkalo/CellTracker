function [seprate_membrane_intensity_all,seprate_membrane]=get_local_membranemasks(nuc_dat1,nuc_dat2,binary_mask,seprate_mask2_clean,binary_mask1,nbrh_pxl_sz,im_bkgd_subtracted,local2)
seprate_membrane_intensity_all=struct;
seprate_membrane_all=struct;

for ii=1:size(binary_mask,2)
cellID1 = (1:size(nuc_dat1(ii).xy,1))';%
cellID2 = (1:size(nuc_dat2(ii).xy,1))';%
celltype1 = cat(2,nuc_dat1(ii).xy, cellID1);%
celltype2 = cat(2,nuc_dat2(ii).xy, cellID2);%
seprate_membrane_intensity=struct;

pixels_membrane = regionprops(binary_mask(ii).dat,'PixelIdxList');% pixels of the area around the cell,ii=image index
pixels_membrane1 = cat(1,pixels_membrane.PixelIdxList);
% now get the intersecion of the pixel area around the given cell with the
% membrane pixels in that image (ii)
%close all
seprate_membrane=struct;
combine_masks(ii).binary = seprate_mask2_clean(ii).dat + binary_mask1(ii).dat;

for cellID=1:size(celltype2,1)
    %disp(cellID)
T = zeros(size(binary_mask(ii).dat));    
T(round(celltype2(cellID,1)),round(celltype2(cellID,2)))=1;
T_1 = imdilate(T,strel('disk',round(nbrh_pxl_sz*0.6)));% such that only look at immediate membranes
T_2 = T_1';
%stats_dilarea = regionprops(T_2,'PixelIdxList');% pixels of the area around the cell
newImg = (size(binary_mask(ii).dat));
newImg = binary_mask(ii).dat & T_2; % pixels of the membrane mask that are only in the neighborhood of a cell
stats = regionprops(newImg,im_bkgd_subtracted(ii).im,'MeanIntensity');
seprate_membrane_intensity(cellID).int = mean(cat(1,stats.MeanIntensity));
seprate_membrane(cellID).maks = newImg;
if mod(cellID,80) == 0
figure(2),imshowpair(combine_masks(ii).binary,newImg);hold on
figure(2),plot(celltype2(cellID,1),celltype2(cellID,2),'pb','MarkerFaceColor','b'); hold on
if ~isempty(local2(ii).dat(cellID).same)
figure(2),plot(local2(ii).dat(cellID).same(:,1),local2(ii).dat(cellID).same(:,2),'*b'); hold on
title('Blue filled = current cell;Neighbors:  Blue = same cell type; red = unlike cells')
end
if ~isempty(local2(ii).dat(cellID).other)
figure(2),plot(local2(ii).dat(cellID).other(:,1),local2(ii).dat(cellID).other(:,2),'*r'); hold on
title('Blue filled = current cell;Neighbors:  Blue = same cell type; red = unlike cells')
end
end
end
seprate_membrane_intensity_all(ii).dat=seprate_membrane_intensity;
seprate_membrane_all(ii).dat=seprate_membrane;
end
disp('done');
end