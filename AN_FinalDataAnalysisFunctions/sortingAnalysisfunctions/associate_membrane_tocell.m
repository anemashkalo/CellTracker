% roughly asssociate the membrane of the cell to a cell
% here assumed that cell types are labelled with DAPI (both cell types have dapi) and CFP (only
% one cell type has this label). 
% Then use the full membrane masks obtained from 'analyze_membrane_stain.m'
% for separation of membrane pixels (see codes below)
maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-02-CAMdyn_Kong_EcadCY5_NcadRFP_CPFundiff_DAPI_RFP_CY5_CFP\maxProjections';
str_type1='_CFP_Probabilities.h5';
str_type2='_DAPI_Probabilities.h5';

[mask1_out_h5,mask2_out_h5,mask_out1,mask_out2,all_h5_type1,all_h5_type2]=get_masks_twocelltypes(maxproj_dir,str_type1,str_type2);


%% threshold the masks and separate the masks for each cell type
prob_thresh = 0.8;
small_stuff = 80;
small_stuff2 = 200;
[binary_mask1,seprate_mask2_clean]=get_celltype_masks(mask_out1,mask_out2,prob_thresh,small_stuff,small_stuff2);
% test output:

% for jj=1:size(binary_mask1,2)
% figure, imshowpair(seprate_mask2_clean(jj).dat,binary_mask1(jj).dat);
% end

%% get the centroids for the cell type specific masks, also get the membrane stats (pixelIdLists)
membrane_matfile = 'MembraneEcad_quantification.mat';
close all
load(membrane_matfile)
[newstr,raw_im,nuc_dat1,nuc_dat2]=get_celltypes_centroids(all_h5_type2,mask_out2,binary_mask1,seprate_mask2_clean);
% test output:
% for k=1:size(binary_mask,2)
% figure,imshowpair(seprate_mask2_clean(k).dat,binary_mask(k).dat);hold on
% plot(nuc_dat1(k).xy(:,1),nuc_dat1(k).xy(:,2),'pr'); hold on
% plot(nuc_dat2(k).xy(:,1),nuc_dat2(k).xy(:,2),'*b'); hold on
% end
membrane_masks= struct;
membrane_masks = binary_mask;
%save('celltype_centroids_masks','nuc_dat1','nuc_dat2','prob_thresh','seprate_mask2_clean','binary_mask1','membrane_masks')
%% now find the local neighborhood of each cell and tell whether these are like or unlike cell
% type relative to the current cell

nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)
nbrh_pxl_sz = 50;%50
% nuc_dat1 = cfp cells
% nuc_dat2 = non-cfp cells

[local2]= get_local_allimges(binary_mask,nuc_dat1,nuc_dat2,nbrh_pxl_sz,seprate_mask2_clean,binary_mask1);

%save('celltype_centroids_masks','local2','nbrh_pxl_sz','-append');% rename the


%% now get the intersection between pixels in the membrane mask with the
% area of pixels surroinding the current cell generated as the dilation of
% the current nuclear centroid
nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)

[seprate_membrane_intensity_all,seprate_membrane_all]=get_local_membranemasks(nuc_dat1,nuc_dat2,binary_mask,seprate_mask2_clean,binary_mask1,nbrh_pxl_sz,im_bkgd_subtracted,local2);
%save('celltype_centroids_masks','seprate_membrane_intensity_all','seprate_membrane_all','-append');% 

%% now split the membrane surrounding the diven cell into regions being the interfaces between neighhboring cells
% using the same approach as used above for obtaining the membrane in the
% neighborhood of the cell

nuc_data_matfle = 'celltype_centroids_masks.mat';
load(nuc_data_matfle);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)

% nuc_dat1 = cfp cells
% nuc_dat2 = non-cfp cells
scale_nbrh = 0.55;
round(nbrh_pxl_sz*scale_nbrh)
str_strel = 'disk';%disk

[membrane_interface_all,membrane_assign_movie]=...
    assign_membrane_segment_to_interface_type(nbrh_pxl_sz,scale_nbrh,str_strel,nuc_dat1,nuc_dat2,seprate_mask2_clean,binary_mask1,newstr,seprate_membrane_all,local2,binary_mask,im_bkgd_subtracted);

% save('celltype_centroids_masks','membrane_interface_all','scale_nbrh','-append');
%%
close all
h = figure(1);h.Colormap = jet; %
movie(h,membrane_assign_movie,1,3);%  last argment: frame rate
% save the movie
tmp_mov = VideoWriter(['12hr_tpnt_MembraneInterfaces_5cells.avi']);
tmp_mov.FrameRate = 1;
% open(tmp_mov);
% writeVideo(tmp_mov,membrane_assign_movie);
% close(tmp_mov);

%% plot eCad vs membrane interface type

membramne_data_file = 'celltype_centroids_masks.mat';
load(membramne_data_file);
membrane_matfile = 'MembraneEcad_quantification.mat';
load(membrane_matfile)
str_full=cell(1,size(newstr,2));
close all
like_interface=struct;
unlike_interface=struct;
for ii=1:size(newstr,2)
str_full{ii} = newstr(ii).imgname(3:end-10);

like_interface(ii).pxlintentsity = cat(1,membrane_interface_all(ii).dat.homotypic);
disp(['Mean in like interfaces: ' num2str(mean(like_interface(ii).pxlintentsity(~isnan(like_interface(ii).pxlintentsity))))])
unlike_interface(ii).pxlintentsity = cat(1,membrane_interface_all(ii).dat.heterotypic);
disp(['Mean in Un_like interfaces: ' num2str(mean(unlike_interface(ii).pxlintentsity(~isnan(unlike_interface(ii).pxlintentsity))))])

figure(ii),histogram(like_interface(ii).pxlintentsity,'binwidth',10,'normalization','probability');hold on
figure(ii),histogram(unlike_interface(ii).pxlintentsity,'binwidth',10,'normalization','probability');box on
legend(['Like interfaces; Mean =: ' num2str(mean(like_interface(ii).pxlintentsity(~isnan(like_interface(ii).pxlintentsity))))],['Unlike interfaces; Mean =: ' num2str(mean(unlike_interface(ii).pxlintentsity(~isnan(unlike_interface(ii).pxlintentsity))))]);
ylim([0 1])
xlabel('E-cadherin mean intensity, a.u.');
ylabel('Freguency');
ylim([0 0.35])
xlim([0 200])
h = figure(ii);
%h.CurrentAxes.XTick = 1:size(newstr,2);
%h.CurrentAxes.XTickLabel=str_full;
h.CurrentAxes.LineWidth = 1.5;
h.CurrentAxes.FontSize = 12;
box on
figure(ii),title(['Ecad intensity in section of the membrane interface  at ' str_full{ii}] );

end


% figure(2),plot(like_interface,'rp'); hold on
% plot(unlike_interface,'bp'); hold on
% legend('Like interfaces','Unlike interfaces');




