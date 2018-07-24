% roughly asssociate the membrane of the cell to a cell
% here assumed that cell types are labelled with DAPI (both cell types have dapi) and CFP (only
% one cell type has this label). Code below separates the nuclear masks into two separate
% channles such that each cell type has it's own mask
% Then use the full membrane masks obtained from 'analyze_membrane_stain.m'
% for separation of membrane pixels (see comments below)
maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-02-CAMdyn_Kong_EcadCY5_NcadRFP_CPFundiff_DAPI_RFP_CY5_CFP\maxProjections';
ff = dir(maxproj_dir);
all_h5_type1 = struct;
all_h5_type2 = struct;
nuc_mask1 = struct;% for cell type 1
nuc_mask2 = struct;% for cell type 2
mask1_out_h5 = struct;
mask_out1 = struct;
mask2_out_h5 = struct;
mask_out2 = struct;
q = 1;
s = 1;
for ii =1:size(ff,1)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'_CFP_Probabilities.h5'))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5_type1(q).file = ff(ii).name; 
nuc_mask1(q).raw =all_h5_type1(q).file;
mask1_out_h5(q).raw = h5read(nuc_mask1(q).raw,'/exported_data');
nuc_lbl = 2;
mask_out1(q).mask = squeeze(mask1_out_h5(q).raw(nuc_lbl,:,:));
q = q+1;
    end
if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'_DAPI_Probabilities.h5'))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5_type2(s).file = ff(ii).name; 
nuc_mask2(s).raw =all_h5_type2(s).file;
mask2_out_h5(s).raw = h5read(nuc_mask2(s).raw,'/exported_data');
nuc_lbl = 2;
mask_out2(s).mask = squeeze(mask2_out_h5(s).raw(nuc_lbl,:,:));
s = s+1;
    end
end

%% threshold the masks and separate the masks for each cell type
prob_thresh = 0.8;
binary_mask1 = struct;
binary_mask2 = struct;
small_stuff = 80;
for jj=1:size(mask_out1,2)
%close all
binary_mask1(jj).dat = mask_out1(jj).mask >prob_thresh;
% filter area
binary_mask1(jj).dat = bwareaopen(binary_mask1(jj).dat,small_stuff)';% ilastik returns transposed image, need to T back
figure,imshowpair(mask_out1(jj).mask(:,:)',binary_mask1(jj).dat);
%figure, imshow(binary_mask(jj).dat,[]);
end
for jj=1:size(mask_out2,2)
%close all
binary_mask2(jj).dat = mask_out2(jj).mask >prob_thresh;
% filter area
binary_mask2(jj).dat = bwareaopen(binary_mask2(jj).dat,small_stuff)';% ilastik returns transposed image, need to T back
figure,imshowpair(mask_out2(jj).mask(:,:)',binary_mask2(jj).dat);
%figure, imshow(binary_mask(jj).dat,[]);
end
% separate the cell types; mask2 corresponds to mask in dapi (all cells are
% present) while mask1 is one seprate cell type1; this will make mask2
% only containing cell type2  
small_stuff2 = 200;
seprate_mask2 = struct;
seprate_mask2_clean = struct;
I_tmp = [];
for jj=1:size(mask_out2,2)
figure,imshowpair(binary_mask2(jj).dat,binary_mask1(jj).dat);
I_tmp = binary_mask2(jj).dat;
seprate_mask2(jj).dat = I_tmp&~binary_mask1(jj).dat;% these are seprately cell type 2muc mask
figure, imshow(seprate_mask2(jj).dat,[]);
% need to clean up the small stuff again (left after removing pixels related to cell type 1)
seprate_mask2_clean(jj).dat = bwareaopen(seprate_mask2(jj).dat,small_stuff2);% ilastik returns transposed image, need to T back


figure, imshowpair(seprate_mask2_clean(jj).dat,binary_mask1(jj).dat);
end

%% get the centroids for the cell type specific masks, also get the membrane stats (pixelIdLists)
membrane_matfile = 'MembraneEcad_quantification.mat';
close all
load(membrane_matfile)
ff = dir(maxproj_dir);
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

figure,imshowpair(seprate_mask2_clean(ii).dat,binary_mask(ii).dat);hold on
plot(nuc_dat1(ii).xy(:,1),nuc_dat1(ii).xy(:,2),'pr'); hold on
plot(nuc_dat2(ii).xy(:,1),nuc_dat2(ii).xy(:,2),'*b'); hold on

end

 %save('celltype_centroids_masks','nuc_dat1','nuc_dat2','prob_thresh','seprate_mask2_clean','binary_mask1')
%% now find the local neighborhood of each cell and tell whether these are like or unlike cell
% type relative to the current cell
%ii will be the loop over images later (experimental conditions)
% clear celltype1 
% clear celltype2

nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)
ii = 1;
nbrh_pxl_sz = 50;%50
combine_masks = struct;
local=struct;
% flip the assigment of celltype1 and celltype2, in order to look at the
% other cell type without chaning the code.
% nuc_dat1 = cfp cells
% nuc_dat2 = non-cfp cells
cellID1 = (1:size(nuc_dat1(ii).xy,1))';%
cellID2 = (1:size(nuc_dat2(ii).xy,1))';%
celltype1 = cat(2,nuc_dat1(ii).xy, cellID1);%
celltype2 = cat(2,nuc_dat2(ii).xy, cellID2);%

for cellID=1:size(celltype2,1)% will need to loop over cells
[like_cells,unlike_cells,no_neighbors] = find_neighbors_simulation(cellID,celltype2,celltype1,nbrh_pxl_sz);
% the cellID is from the cell type that follows this argument (cells from
% celltype1) in this case
disp(cellID)
if mod(cellID,20)==0
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

%save('celltype_centroids_masks','local','ii','nbrh_pxl_sz','-append');% rename the

%localstructure specific to the image analysed(ii)

%% now get the intersection between pixels in the membrane mask with the
% area of pixels surroinding the current cell generated as the dilation of
% the current nuclear centroid
nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)
% flip the assigment of celltype1 and celltype2, in order to look at the
% other cell type without chaning the code.
% nuc_dat1 = cfp cells
% nuc_dat2 = non-cfp cells
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
    disp(cellID)
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
if mod(cellID,20) == 0
figure(2),imshowpair(combine_masks(ii).binary,newImg);hold on
figure(2),plot(celltype2(cellID,1),celltype2(cellID,2),'pb','MarkerFaceColor','b'); hold on
if ~isempty(local(cellID).same)
figure(2),plot(local(cellID).same(:,1),local(cellID).same(:,2),'*b'); hold on
title('Blue filled = current cell;Neighbors:  Blue = same cell type; red = unlike cells')
end
if ~isempty(local(cellID).other)
figure(2),plot(local(cellID).other(:,1),local(cellID).other(:,2),'*r'); hold on
title('Blue filled = current cell;Neighbors:  Blue = same cell type; red = unlike cells')
end
end
end
disp('done');
%save('celltype_centroids_masks','seprate_membrane','seprate_membrane_intensity','-append');% rename the
% same make the membrane intensity masks specific to image (ii)
%% now split the membrane surrounding the diven cell into regions being the interfaces between neighhboring cells
% using the same approach as used above for obtaining the membrane in the
% neighborhood of the cell

nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)
% flip the assigment of celltype1 and celltype2, in order to look at the
% other cell type without chaning the code.
% nuc_dat1 = cfp cells
% nuc_dat2 = non-cfp cells
cellID1 = (1:size(nuc_dat1(ii).xy,1))';%
cellID2 = (1:size(nuc_dat2(ii).xy,1))';%
celltype1 = cat(2,nuc_dat1(ii).xy, cellID1);%
celltype2 = cat(2,nuc_dat2(ii).xy, cellID2);%
combine_masks(ii).binary = seprate_mask2_clean(ii).dat + binary_mask1(ii).dat;

seprate_membrane_intensity = struct;
membrane_interface = struct;
q = 1;
s = 1;
qq = 1;
membrane_assign_movie = struct('cdata',[],'colormap',[]);
close all

newStr = replace(newstr(ii).imgname,'CY5','CY5');
img = imread(newStr);
img = imadjust(img,[0 0.4]);%stretchlim(img)
img = im2uint8(img);
newStr2 = replace(newstr(ii).imgname,'CY5','DAPI');
img2 = imread(newStr2);
img2 = imadjust(img2,[0 0.6]);%stretchlim(img2)
img3 = max(img,img2);
scale_nbrh = 0.55;
round(nbrh_pxl_sz*scale_nbrh)
str_strel = 'disk';

for cellID=1:size(celltype2,1)
    disp(['left to process  ' num2str(size(celltype2,1)-cellID)])    
% nowloop over neighbors and get the membrane pixels corresponding to interace
% between neighbor then 
separate_membrane_mask_tmp = seprate_membrane(cellID).maks;
if ~isempty(local(cellID).same)
for jj=1:size(local(cellID).same,1)
%figure(4), imshow(separate_membrane_mask_tmp) ;hold on;% get rid of the 
    
T = zeros(size(binary_mask(ii).dat)); 
T(round(local(cellID).same(jj,1)),round(local(cellID).same(jj,2)))=1;
T_1 = imdilate(T,strel(str_strel,round(nbrh_pxl_sz*scale_nbrh)));% such that only look at immediate membranes
T_2 = T_1';
newImg_tmp = zeros(size(binary_mask(ii).dat));
newImg_tmp = separate_membrane_mask_tmp & T_2; % pixels of the membrane mask that are only in the neighborhood of a given cell
% and its local neighbor
separate_membrane_mask_tmp = separate_membrane_mask_tmp & ~newImg_tmp;
% assigned membrane segment to avoid contribution of same pixels to
% different 
if sum(sum(newImg_tmp)) >0
stats = regionprops(newImg_tmp,im_bkgd_subtracted(ii).im,'MeanIntensity');
seprate_membrane_intensity(cellID).int = mean(cat(1,stats.MeanIntensity));
membrane_interface(q).homotypic = seprate_membrane_intensity(cellID).int;
q=q+1;
end
% look at the membrane that was selected
if mod(cellID,20) == 0   
    %,'ColorChannels','red-cyan'
figure(3), imshowpair(newImg_tmp,img,'method','blend');hold on%seprate_membrane(cellID).maks
plot(celltype2(:,1),celltype2(:,2),'g.','MarkerSize',7); hold on
plot(celltype1(:,1),celltype1(:,2),'r.','MarkerSize',7); hold on
plot(celltype2(cellID,1),celltype2(cellID,2),'pw','MarkerFaceColor','b','MarkerSize',5); hold on
title('Blue filled = current cell;Neighbors:  green = same cell type; red = unlike cells')
if ~isempty(local(cellID).same)
figure(3),plot(local(cellID).same(jj,1),local(cellID).same(jj,2),'.g','MarkerSize',15); hold on
title('Blue filled = current cell;Neighbors: green = same cell type; red = unlike cells')
end
 h3 = figure(3);
 %h3.Position
        %h3.Position =[0 0 1024 1024];%%560 420 this is the size that the function 'movie' supports
        membrane_assign_movie(qq) = getframe(gcf,[0 0 686   605]);%h3,h3.Position getframe(h3)
       % h3.Units = 'normalized';
        qq = qq+1;
end
end
end
% same for the heretotypic interface
if ~isempty(local(cellID).other)
for jj=1:size(local(cellID).other,1)
%figure(4), imshow(separate_membrane_mask_tmp) ;hold on;% get rid of the 

T = zeros(size(binary_mask(ii).dat)); 
T(round(local(cellID).other(jj,1)),round(local(cellID).other(jj,2)))=1;
T_1 = imdilate(T,strel(str_strel,round(nbrh_pxl_sz*scale_nbrh)));% such that only look at immediate membranes
T_2 = T_1';
newImg_tmp = zeros(size(binary_mask(ii).dat));
newImg_tmp = separate_membrane_mask_tmp & T_2; % pixels of the membrane mask that are only in the neighborhood of a given cell
% and its local neighbor
separate_membrane_mask_tmp = separate_membrane_mask_tmp & ~newImg_tmp;
% assigned membrane segment to avoid contribution of same pixels to
% different
if sum(sum(newImg_tmp)) >0
stats = regionprops(newImg_tmp,im_bkgd_subtracted(ii).im,'MeanIntensity');
seprate_membrane_intensity(cellID).int = mean(cat(1,stats.MeanIntensity));
membrane_interface(s).heterotypic = seprate_membrane_intensity(cellID).int; 
s=s+1;    
end
% look at the membrane that was selected
if mod(cellID,20) == 0    
figure(3), imshowpair(newImg_tmp,img,'method','blend');hold on%seprate_membrane(cellID).maks
plot(celltype2(:,1),celltype2(:,2),'g.','MarkerSize',7); hold on
plot(celltype1(:,1),celltype1(:,2),'r.','MarkerSize',7); hold on

plot(celltype2(cellID,1),celltype2(cellID,2),'pw','MarkerFaceColor','b','MarkerSize',5); hold on
title('Blue filled = current cell;Neighbors:  green = same cell type; red = unlike cells')

if ~isempty(local(cellID).other)
figure(3),plot(local(cellID).other(jj,1),local(cellID).other(jj,2),'.r','MarkerSize',15); hold on
title('Blue filled = current cell;Neighbors:  green = same cell type; red = unlike cells')
end

h3 = figure(3);        
        %h3.Position =[0 0 1024 1024];%%560 420 this is the size that the function 'movie' supports
        membrane_assign_movie(qq) = getframe(gcf,[0 0 686   605]);%h3,h3.Position getframe(h3)
        %h3.Units = 'normalized';
        qq = qq+1;
end
end
end
end
% save('celltype_centroids_masks','membrane_interface','scale_nbrh','-append');
% make labelled by the image number ii
%%
close all
h = figure(1);h.Colormap = jet; %
movie(h,membrane_assign_movie,1,3);%  last argment: frame rate
%% save the movie

tmp_mov = VideoWriter(['2hr_time_pnt_MembraneInterfaces.avi']);
tmp_mov.FrameRate = 1;
% open(tmp_mov);
% writeVideo(tmp_mov,membrane_assign_movie);
% close(tmp_mov);



%% plot eCad vs membrane interface type
close all
like_interface = cat(1,membrane_interface.homotypic);
disp(['Mean in like interfaces: ' num2str(mean(like_interface(~isnan(like_interface))))])
unlike_interface = cat(1,membrane_interface.heterotypic);
disp(['Mean in Un_like interfaces: ' num2str(mean(unlike_interface(~isnan(unlike_interface))))])

figure(1),histogram(like_interface,'binwidth',10,'normalization','probability');hold on
histogram(unlike_interface,'binwidth',10,'normalization','probability');box on
legend(['Like interfaces; Mean =: ' num2str(mean(like_interface(~isnan(like_interface))))],['Unlike interfaces; Mean =: ' num2str(mean(unlike_interface(~isnan(unlike_interface))))]);
ylim([0 1])
if ii == 4
title('Ecad intensity in section of the membrane interface, time point: 24 hr' );
end
if ii == 3
title('Ecad intensity in section of the membrane interface, time point: 12 hr' );
end
if ii == 2
title('Ecad intensity in section of the membrane interface, time point: 6 hr' );
end
if ii == 1
title('Ecad intensity in section of the membrane interface, time point: 2 hr' );
end
xlabel('E-cadherin mean intensity, a.u.');
ylabel('Freguency');
ylim([0 0.35])
xlim([0 200])

figure(2),plot(like_interface,'rp'); hold on
plot(unlike_interface,'bp'); hold on
legend('Like interfaces','Unlike interfaces');


%% look at the fractions of unlike cells vs like cells and get the
% meanFluorescence for eCad based on that

loc_fractions_same = zeros(size(celltype2,1),1);%size(celltype2,1),1
loc_fractions_other = zeros(size(celltype2,1),1);%size(celltype2,1),1
loc_intensities = zeros(size(celltype2,1),1);%size(celltype2,1),1
total = zeros(size(celltype2,1),1);%size(celltype2,1),1

for cellID=1:size(celltype2,1)
total(cellID) = size(local(cellID).same,1)+size(local(cellID).other,1);    
loc_fractions_same(cellID) = size(local(cellID).same,1)/total(cellID);
loc_fractions_other(cellID) = size(local(cellID).other,1)/total(cellID);

end
fractions_same=loc_fractions_same(~isnan(loc_fractions_same));
fractions_other=loc_fractions_other(~isnan(loc_fractions_other));
fin_intensity_ecad = cat(1,seprate_membrane_intensity.int);
fin_intensity_ecad = fin_intensity_ecad(~isnan(loc_fractions_other));
figure(1), plot(fin_intensity_ecad,fractions_same,'p');box on;hold on
corrcoef(fin_intensity_ecad(~isnan(fin_intensity_ecad)),fractions_same(~isnan(fin_intensity_ecad)))

xlabel('Ecadherin intensity in the local neighborhood of given cell');
ylabel('Fraction of neighborhood that is the same cells type');
figure(1), plot(fin_intensity_ecad,fractions_other,'p');box on
xlabel('Ecadherin intensity in the local neighborhood of given cell');
ylabel('Fraction of neighborhood that is the unlike cells type');
corrcoef(fin_intensity_ecad(~isnan(fin_intensity_ecad)),fractions_other(~isnan(fin_intensity_ecad)))

total = total(~isnan(loc_fractions_other));
figure(3), plot(fin_intensity_ecad,total,'p');box on;hold on
xlabel('Ecadherin intensity in the local neighborhood of given cell');
ylabel('Total number of neighbors');
corrcoef(fin_intensity_ecad(~isnan(fin_intensity_ecad)),total(~isnan(fin_intensity_ecad)))



