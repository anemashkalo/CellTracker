%% process membrane fluorescnece from singleimages and ilastic masks
maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-31-EcadReporterYFPwithbetacatRFP_liveimaging\MaxProjections';
%maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-04-18-EcadNcad_representative\max_proj_bychannel_betaCatcells_CY5ecad';
ff = dir(maxproj_dir);
all_h5 = struct;
membrane_mask = struct;
mask_out_h5 = struct;
mask_out = struct;
q = 1;
for ii =1:size(ff,1)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'_CY5_Probabilities.h5'))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5(q).file = ff(ii).name; 
membrane_mask(q).raw =all_h5(q).file;
mask_out_h5(q).raw = h5read(membrane_mask(q).raw,'/exported_data');
membr_lbl = 2;
mask_out(q).mask = squeeze(mask_out_h5(q).raw(membr_lbl,:,:));
q = q+1;
end
end

%%
prob_thresh = 0.85;
binary_mask = struct;
small_stuff = 60;
for jj=1:size(mask_out,2)
%close all
binary_mask(jj).dat = mask_out(jj).mask >prob_thresh;
% filter area
binary_mask(jj).dat = bwareaopen(binary_mask(jj).dat,small_stuff)';% ilastik returns transposed image, need to T back
%figure,imshowpair(mask_out(jj).mask(:,:)',binary_mask(jj).dat);
%figure,imshow(binary_mask(jj).dat,[]);
end
%% subtract background and quntify fluorescence
close all
ff = dir(maxproj_dir);
raw_im = struct;
im_bkgd_subtracted = struct;
dat = struct;
newstr = struct;
mean_ecad = zeros(size(all_h5,2),1);
err_ecad= zeros(size(all_h5,2),1);
for ii =1:size(all_h5,2) 
newstr(ii).imgname = [all_h5(ii).file(1:end-17) '.tif' ];
raw_im(ii).dat = imread(newstr(ii).imgname);%im2uint8
[im_bkgd_subtracted(ii).im] = simplebg([],binary_mask(ii).dat,raw_im(ii).dat);
%figure,imshowpair(im_bkgd_subtracted(ii).im,binary_mask(ii).dat);
figure,imshowpair(raw_im(ii).dat,binary_mask(ii).dat);
% get mean fluorescence intensity
stats = regionprops(binary_mask(ii).dat,im_bkgd_subtracted(ii).im,'MeanIntensity','Area','PixelIdxList');
dat(ii).ecadlevels = cat(1,stats.MeanIntensity);
dat(ii).area = cat(1,stats.Area);
%dat(ii).pxls = stats.PixelIdxList;
dat(ii).objects = size(stats,1);
%disp(dat(ii).objects);
mean_ecad(ii,1) = mean(dat(ii).ecadlevels);%./dat(ii).area
err_ecad(ii,1) = std(dat(ii).ecadlevels);
end
prob_thresh_membr = prob_thresh;
%save('MembraneEcad_quantification_BetaCatMembrane','dat','newstr','binary_mask','prob_thresh_membr','im_bkgd_subtracted');
%%
close all
str = cell(1,size(all_h5,2));
%str ={'betaCat prediff only','betaCat pluri only','betaCat pluri, esi diff','betaCat diff, esi pluri'};
for ii =1:size(all_h5,2) 
figure(ii+1), histogram(dat(ii).ecadlevels,'binwidth',500,'normalization','probability'); hold on
legend(str{ii});box on
end
figure(1), plot(mean_ecad,'kp','MarkerSize',12,'MarkerFaceColor','m'); box on
h = figure(1);
h.CurrentAxes.XLim = [0 size(all_h5,2)+1];
h.CurrentAxes.YLim = [0 max(mean_ecad)+500];%
h.CurrentAxes.XTick = [1:size(all_h5,2)];
h.CurrentAxes.LineWidth = 2;
h.CurrentAxes.FontSize = 12;
% h.CurrentAxes.XTickLabel={'betaCat prediff only','betaCat pluri only','betaCat pluri, esi diff','betaCat diff, esi pluri'};
h.CurrentAxes.XTickLabelRotation = -25;
%title('Fixed after 30 hrs of sorting')
%ylabel('Membrane Beta-catenin in reporter cells, a.u.')
% xlabel('Time after mixing the two cell types or unmixed fixation time')



