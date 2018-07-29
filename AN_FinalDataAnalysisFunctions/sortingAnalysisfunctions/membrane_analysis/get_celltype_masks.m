function [binary_mask1,seprate_mask2_clean]=get_celltype_masks(mask_out1,mask_out2,prob_thresh,small_stuff,small_stuff2)
% this function separates the masks for two cell types using the images
%where all cells are DAPI+ and cell type1 is other channel+ (cfp+)
% this function makesmasks exclusively for each cell type(no dapi overlap for both cell types)

binary_mask1 = struct;
binary_mask2 = struct;

for jj=1:size(mask_out1,2)
%close all
binary_mask1(jj).dat = mask_out1(jj).mask >prob_thresh;
% filter area
binary_mask1(jj).dat = bwareaopen(binary_mask1(jj).dat,small_stuff)';% ilastik returns transposed image, need to T back
%figure,imshowpair(mask_out1(jj).mask(:,:)',binary_mask1(jj).dat);
%figure, imshow(binary_mask(jj).dat,[]);
end
for jj=1:size(mask_out2,2)
%close all
binary_mask2(jj).dat = mask_out2(jj).mask >prob_thresh;
% filter area
binary_mask2(jj).dat = bwareaopen(binary_mask2(jj).dat,small_stuff)';% ilastik returns transposed image, need to T back
%figure,imshowpair(mask_out2(jj).mask(:,:)',binary_mask2(jj).dat);
%figure, imshow(binary_mask(jj).dat,[]);
end
% separate the cell types; mask2 corresponds to mask in dapi (all cells are
% present) while mask1 is one seprate cell type1; this will make mask2
% only containing cell type2  

seprate_mask2 = struct;
seprate_mask2_clean = struct;
I_tmp = [];
for jj=1:size(mask_out2,2)
%figure,imshowpair(binary_mask2(jj).dat,binary_mask1(jj).dat);
I_tmp = binary_mask2(jj).dat;
seprate_mask2(jj).dat = I_tmp&~binary_mask1(jj).dat;% these are seprately cell type 2muc mask
%figure, imshow(seprate_mask2(jj).dat,[]);
% need to clean up the small stuff again (left after removing pixels related to cell type 1)
seprate_mask2_clean(jj).dat = bwareaopen(seprate_mask2(jj).dat,small_stuff2);% ilastik returns transposed image, need to T back
%figure, imshowpair(seprate_mask2_clean(jj).dat,binary_mask1(jj).dat);
end
end