%% determine the background images(for all chnnels) for the dataset to run

%run once before checking parameters using the section below

% imagedirectoryname has to be a string 'imagedirectoryname'. need to be
% one dir up from raw images

ff=readMMdirectory(imagedirectoryname);   % 
dims = [ max(ff.pos_x)+1 max(ff.pos_y)+1];
wavenames=ff.chan;
maxims= dims(1)*dims(2);
%generate background image for each channel
    for ii=1:length(wavenames) % get the background image for al channels  
        [minI, meanI]=mkBackgroundImageMM(ff,ii,min(500,maxims));
        bIms{ii}=uint16(2^16*minI);
        nIms{ii}=ones(size(bIms{ii}));
    end
%% script to optimize the segmentation parameters. Can look at a chosen image
% and adjust the parameters in the paramfile, if segmentation is not good. N is a linear index, image number
% need to be one directory up from the actual images folder ( since using
% the readMM2irectory function here)% 
% the parameter values below are mostly what you will need to adjust to
% obtain resonable segmentation. These need to be changed in the paramfile,
% not here

% userParam.nucIntensityLoc     (line 71 in paramfile): this parameter
% strongly depends on the mean fluorescence intensity of the background (use improfile to estimate this value and use it as a starting 
% test value for this parameter and see how the segmentation improves as you change it)

% userParam.gaussRadius         (line 37 in paramfile)
% userParam.gaussSigma          (line 38 in paramfile)

paramfile=('/Users/warmflashlab/CellTracker/paramFiles/setUserParamAN20X_uCOL.m');
close all
N =1;%  
flag = 1;% plot the resulting segmentation on images
[data1,mask1,toshow1]=ANrunOneMM(imagedirectoryname,N,bIms,nIms,'setUserParamAN20X_uCOL','DAPI',flag);%setUserParamAN20X  setUserParamAN20X_uCOL
imcontrast
%% once happy with the parameter set and the segmentationquality run this to process all the images in the 
% directory
% make a more descriptive name for the outputfile name (I usually put the stains in order into the filename, so I remember after a while what was stained for)
runFullTileMM(imagedirectoryname,'outputfile1.mat','setUserParamAN20X_uCOL',1);
runFullTileMM(imagedirectoryname,'outputfile2.mat','setUserParamAN20X_uCOL',1);

disp('done');



