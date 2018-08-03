% analze membrane stain SD confocal movies (supply max projections and ilastik masks)

maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-31-EcadReporterYFPwithbetacatRFP_liveimaging\MaxProjections';
%%
ff = readAndorDirectory(maxproj_dir);
fnm = struct;
small_stuff=300;
prob_thresh=0.55;
delta_t = 30; % minutes
ilastik = struct;
chan_membrane = 1;

chan = ff.w(chan_membrane); %
processed_movies = struct;
for ii =1:size(ff.p,2)
    pos = ff.p(ii);
    membranedat = struct;
   raw_im = struct;
   im_bkgd_subtracted = struct;
   reader = [];
 mask2=[];
mask3=[];

   fnm(ii).name = getAndorFileName(ff,pos,[],[],ff.w(chan_membrane));      %getAndorFileName(files,pos,time,z,w)
   ilastik(ii).name = [fnm(ii).name(1:end-4) '_Probabilities.h5'];
   disp(['loading: ' num2str(fnm(ii).name)])        
reader = bfGetReader(fnm(ii).name);
nz=reader.getSizeZ;
nT = reader.getSizeT;
% get masks
mask2 = readIlastikProbMask(ilastik(ii).name,prob_thresh);
mask3 = bwareaopen(mask2,small_stuff);
% now subtract background from images and apply masks to get
% MeanPixelIntensity in membrane
stats = [];
for time=1:nT
iPlane=reader.getIndex(nz - 1, chan_membrane -1, time - 1) + 1;
        raw_im(time).dat=bfGetPlane(reader,iPlane);        
        im_bkgd_subtracted(time).im = simplebg([],mask3(:,:,time),raw_im(time).dat); 
       % figure(1),imshowpair(im_bkgd_subtracted(time).im,mask3(:,:,time));
        stats = regionprops(mask3(:,:,time),im_bkgd_subtracted(time).im,'MeanIntensity','Area');
membranedat(time).meanpxlint = mean(cat(1,stats.MeanIntensity));
membranedat(time).err = std(cat(1,stats.MeanIntensity));
membranedat(time).objects = size(stats,1);
end
processed_movies(ii).img = im_bkgd_subtracted;
processed_movies(ii).pxlint = cat(1,membranedat.meanpxlint);
processed_movies(ii).err = cat(1,membranedat.err);
end 
%save('Ecad_analysis','processed_movies','ff','delta_t','prob_thresh');
%% plot the obtained membrane stain quantification
matfile = 'Ecad_analysis.mat';
load(matfile);
close all
max_ecad=[];
C = {'b','b','r','r','m','m','c','c'};
for ii=1:size(processed_movies,2)
figure(2),p1(ii)= plot(processed_movies(ii).pxlint,'-p','Color',C{ii}); hold on
box on
h2 = figure(2);
h2.CurrentAxes.LineWidth = 1.5;
h2.CurrentAxes.FontSize=12;
h2.CurrentAxes.XTick = (1:5:size(processed_movies(ii).pxlint,1));
h2.CurrentAxes.XTickLabel = (1:5:size(processed_movies(ii).pxlint,1))*delta_t/60;
max_ecad(ii) = max(processed_movies(ii).pxlint);
end
xlabel('time, hrs')
ylabel('Ecadherin pixel intensity in membrane of cells, a.u.')
title('Cell Sorting, Ecad reporter cell line')
ylim([100 max(max_ecad)]);
legend([p1(1) p1(3) p1(5) p1(7)],{'Pre-differentaited unmixed','Pluripotent unmixed','Prediff, mixed with pluri','Pluri, mixed with prediff'});


