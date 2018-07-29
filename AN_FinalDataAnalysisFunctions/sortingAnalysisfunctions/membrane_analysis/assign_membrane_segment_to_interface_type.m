function [membrane_interface_all,membrane_assign_movie]=assign_membrane_segment_to_interface_type(nbrh_pxl_sz,scale_nbrh,str_strel,nuc_dat1,nuc_dat2,seprate_mask2_clean,binary_mask1,newstr,seprate_membrane_all,local2,binary_mask,im_bkgd_subtracted)
% this function assigns the local cell membrane to segments,coming from
% homo- or heterotypic interface (based on the local cell neighborhood)

membrane_interface_all=struct;
membrane_assign_movie = struct('cdata',[],'colormap',[]);

for ii=1:size(binary_mask1,2)
seprate_membrane_intensity = struct;
membrane_interface = struct;
separate_membrane_mask_tmp=[];
cellID1 = (1:size(nuc_dat1(ii).xy,1))';%
cellID2 = (1:size(nuc_dat2(ii).xy,1))';%
celltype1 = cat(2,nuc_dat1(ii).xy, cellID1);%
celltype2 = cat(2,nuc_dat2(ii).xy, cellID2);%
combine_masks(ii).binary = seprate_mask2_clean(ii).dat + binary_mask1(ii).dat;
close all
q = 1;
s = 1;
qq = 1;
newStr = newstr(ii).imgname;% needs to be membrane image 
img = imread(newStr);
img = imadjust(img,[0 0.7]);%stretchlim(img)
img = im2uint8(img);
% if want to look at DPI channel instead, or any other channel
% newStr2 = replace(newstr(ii).imgname,'CY5','DAPI');
for cellID=1:size(celltype2,1)
    disp(['left to process  ' num2str(size(celltype2,1)-cellID)])    
% now loop over neighbors and get the membrane pixels corresponding to interace
% between neighbor then 
separate_membrane_mask_tmp = seprate_membrane_all(ii).dat(cellID).maks;
if ~isempty(local2(ii).dat(cellID).same)
for jj=1:size(local2(ii).dat(cellID).same,1)
%figure(4), imshow(separate_membrane_mask_tmp) ;hold on;%   
T = zeros(size(binary_mask(ii).dat)); 
T(round(local2(ii).dat(cellID).same(jj,1)),round(local2(ii).dat(cellID).same(jj,2)))=1;
T_1 = imdilate(T,strel(str_strel,round(nbrh_pxl_sz*scale_nbrh)));% such that only look at immediate membranes
T_2 = T_1';
newImg_tmp = zeros(size(binary_mask(ii).dat));
newImg_tmp = separate_membrane_mask_tmp & T_2; % pixels of the membrane mask that are only in the neighborhood of a given cell
% and its local neighbor
separate_membrane_mask_tmp = separate_membrane_mask_tmp & ~newImg_tmp;

if sum(sum(newImg_tmp)) >0
stats = regionprops(newImg_tmp,im_bkgd_subtracted(ii).im,'MeanIntensity');
seprate_membrane_intensity(cellID).int = mean(cat(1,stats.MeanIntensity));
membrane_interface(q).homotypic = seprate_membrane_intensity(cellID).int;
q=q+1;
% look at the membrane that was selected
if mod(cellID,50) == 0   
    %,'ColorChannels','red-cyan'  'blend'
figure(3), imshowpair(img,(newImg_tmp+seprate_membrane_all(ii).dat(cellID).maks),'method','falsecolor');hold on%  img seprate_membrane(cellID).maks
 plot(celltype2(:,1),celltype2(:,2),'c.','MarkerSize',7); hold on
 plot(celltype1(:,1),celltype1(:,2),'r.','MarkerSize',7); hold on
 plot(celltype2(cellID,1),celltype2(cellID,2),'pw','MarkerFaceColor','b','MarkerSize',5); hold on
title('Blue filled = current cell;Neighbors:  cyan = same cell type; red = unlike cells')
if ~isempty(local2(ii).dat(cellID).same)
figure(3),plot(local2(ii).dat(cellID).same(jj,1),local2(ii).dat(cellID).same(jj,2),'.c','MarkerSize',15); hold on
title('Blue filled = current cell;Neighbors: cyan = same cell type; red = unlike cells')
 h3 = figure(3);
        %membrane_assign_movie(qq) = getframe(gcf,h3.Position);%h3,h3.Position getframe(h3)
        membrane_assign_movie(qq) = getframe(gcf,[0 0 686   605]);%h3,h3.Position getframe(h3)       
        % h3.Units = 'normalized';
        qq = qq+1;
end
end
end
end
% same for the heretotypic interface
if ~isempty(local2(ii).dat(cellID).other)
for jj=1:size(local2(ii).dat(cellID).other,1)
T = zeros(size(binary_mask(ii).dat)); 
T(round(local2(ii).dat(cellID).other(jj,1)),round(local2(ii).dat(cellID).other(jj,2)))=1;
T_1 = imdilate(T,strel(str_strel,round(nbrh_pxl_sz*scale_nbrh)));% such that only look at immediate membranes
T_2 = T_1';
newImg_tmp = zeros(size(binary_mask(ii).dat));
newImg_tmp = separate_membrane_mask_tmp & T_2; % pixels of the membrane mask that are only in the neighborhood of a given cell
% and its local neighbor
separate_membrane_mask_tmp = separate_membrane_mask_tmp & ~newImg_tmp;
% remove the assigned membrane segment from the total local membrane mask to avoid contribution of same pixels to
% different interfae types (avoid overlap)
if sum(sum(newImg_tmp)) >0
stats = regionprops(newImg_tmp,im_bkgd_subtracted(ii).im,'MeanIntensity');
seprate_membrane_intensity(cellID).int = mean(cat(1,stats.MeanIntensity));
membrane_interface(s).heterotypic = seprate_membrane_intensity(cellID).int; 
s=s+1;    
% look at the membrane that was selected for given cell 
%newImg_tmp_crop = imcrop(newImg_tmp,[(round(celltype2(cellID,1))-shift(1)) (round(celltype2(cellID,2))-shift(2)) sz_crop(1) sz_crop(2)]);
if mod(cellID,50) == 0    
figure(3), imshowpair(img,(newImg_tmp+seprate_membrane_all(ii).dat(cellID).maks),'method','falsecolor');hold on%seprate_membrane(cellID).maks
plot(celltype2(:,1),celltype2(:,2),'c.','MarkerSize',7); hold on
plot(celltype1(:,1),celltype1(:,2),'r.','MarkerSize',7); hold on
plot(celltype2(cellID,1),celltype2(cellID,2),'pw','MarkerFaceColor','b','MarkerSize',5); hold on
title('Blue filled = current cell;Neighbors:  cyan = same cell type; red = unlike cells')
if ~isempty(local2(ii).dat(cellID).other)
figure(3),plot(local2(ii).dat(cellID).other(jj,1),local2(ii).dat(cellID).other(jj,2),'.r','MarkerSize',15); hold on
title('Blue filled = current cell;Neighbors:  cyan = same cell type; red = unlike cells')
end
h3 = figure(3);        
        %membrane_assign_movie(qq) = getframe(gcf,h3.Position);
        membrane_assign_movie(qq) = getframe(gcf,[0 0 686   605]);%h3,h3.Position getframe(h3)
        %h3.Units = 'normalized';
        qq = qq+1;
end
end
end
end
end
end
membrane_interface_all(ii).dat = membrane_interface;
end


end