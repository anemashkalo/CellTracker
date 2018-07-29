function [local2]= get_local_allimges(binary_mask,nuc_dat1,nuc_dat2,nbrh_pxl_sz,seprate_mask2_clean,binary_mask1)
% getlocal neighborhood of each cell in each image within the 
% distance nbrh_pxl_sz pizels of the given cell
% the given cell is the cell from arguemnt 3 in the input 
combine_masks = struct;
local2 = struct;
for ii = 1:size(binary_mask,2)
local=struct;
cellID1 = (1:size(nuc_dat1(ii).xy,1))';%
cellID2 = (1:size(nuc_dat2(ii).xy,1))';%
celltype1 = cat(2,nuc_dat1(ii).xy, cellID1);%
celltype2 = cat(2,nuc_dat2(ii).xy, cellID2);%
like_cells=[];
unlike_cells=[];
for cellID=1:size(celltype2,1)% will need to loop over cells
[like_cells,unlike_cells,no_neighbors] = find_neighbors_simulation(cellID,celltype2,celltype1,nbrh_pxl_sz);
% the cellID is from the cell type that follows this argument (cells from
% celltype1) in this case
disp(cellID)
if mod(cellID,50)==0
combine_masks(ii).binary = seprate_mask2_clean(ii).dat + binary_mask1(ii).dat;
figure(2),imshowpair( seprate_mask2_clean(ii).dat,binary_mask1(ii).dat);hold on
if ~isempty(like_cells)
figure(2),plot(celltype2(cellID,1),celltype2(cellID,2),'pb','MarkerFaceColor','b'); hold on
figure(2),plot(like_cells(:,1),like_cells(:,2),'*b'); hold on
end
if ~isempty(unlike_cells)
figure(2),plot(unlike_cells(:,1),unlike_cells(:,2),'*r'); hold on
title('Masks for each cell type separately')
end
figure(2),imshowpair(combine_masks(ii).binary,binary_mask(ii).dat);hold on
plot(celltype2(cellID,1),celltype2(cellID,2),'pb','MarkerFaceColor','b'); hold on
% plot(like_cells(:,1),like_cells(:,2),'*b'); hold on
% plot(unlike_cells(:,1),unlike_cells(:,2),'*r'); hold on
title('Blue filled = current cell;Neighbors:  Blue = same cell type; red = unlike cells')
end
local(cellID).same=like_cells;
local(cellID).other=unlike_cells;
end
local2(ii).dat = local;

end


end
