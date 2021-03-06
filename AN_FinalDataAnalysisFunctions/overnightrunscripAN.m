%%
% run the data analysis overnight
% 

% low  bmp pSmad1-Smad4 corr, use the empty eimage taken in each channel as
% the background (instead of calculating)

% runFullTileMM('pluriControl_4dMT','pluriControl_4dMT.mat','setUserParamAN20X',1);
% runFullTileMM('diffControl_2dBmpSb','diffControl_2dBmpSb.mat','setUserParamAN20X',1);
% runFullTileMM('1d_BmpSb_3d_MT','1d_BmpSb_3d_MT.mat','setUserParamAN20X',1);
% runFullTileMM('2d_BmpSb_2d_MT','2d_BmpSb_2d_MT.mat','setUserParamAN20X',1);
runFullTileMM('3d_BmpSb_1d_MT','3d_BmpSb_1d_MT.mat','setUserParamAN20X',1);
runFullTileMM('4d_BmpSb','4d_BmpSb.mat','setUserParamAN20X',1);

disp('done');

  
%% split the ibidi .btfs
flag = 1;
MMdirec1 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\betaCat_Diff_Ecad';
MMdirec2 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\betaCat_Diff_Ncad';
MMdirec3 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\betaCat_Pluri_Ecad';
MMdirec4 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\betaCat_Pluri_Ncad';
MMdirec5 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\esiDiff_betacatPluri_Ecad';
MMdirec6 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\esiDiff_betacatPluri_Ncad';
MMdirec7 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\esiPluri_betacatDiff_Ecad';
MMdirec8 = 'E:\allSortingData\2018-04-18-EcadNcad_8wellDGRCchan\esiPluri_betacatDiff_Ncad';

fn= {'betaCat_Diff_Ecad.btf','betaCat_Diff_Ncad.btf','betaCat_Pluri_Ecad.btf'...
    ,'betaCat_Pluri_Ncad.btf','esiDiff_betacatPluri_Ecad.btf','esiDiff_betacatPluri_Ncad.btf'...
    ,'esiPluri_betacatDiff_Ecad.btf','esiPluri_betacatDiff_Ncad.btf'};
chan = {'DAPI','GFP','RFP','CY5'};% 
mm = {MMdirec1,MMdirec2,MMdirec3,MMdirec4,MMdirec5,MMdirec6,MMdirec7,MMdirec8};
cellAdhmol = cell(1,size(mm,2));
for k=1:size(mm,2)
   cellAdhmol{k} = templateSplitOlympData(mm{k},chan,fn{k},flag);
end
disp('split and saved all');
save('cellAdhmol');

%% run colony grouping only 

outfile ='300ngmlWNT3a_uCol.mat'; 
direc = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/2_NO_QUADRANTS_goodData(esi017Cells)/2017-01-27-WNT3aDosesandActinuCol/300ngmlWNT3a_uCol';
direc2 = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/2_NO_QUADRANTS_goodData(esi017Cells)/2017-01-27-WNT3aDosesandActinuCol/outfiles_group130px';
paramfile=('/Users/warmflashlab/CellTracker/paramFiles/setUserParamAN20X_uCOL.m');
run(paramfile);
ff = readMMdirectory(direc);
load([direc filesep outfile],'bIms','nIms','dims');
    [colonies, peaks]=peaksToColonies([direc filesep outfile]);
    plate1=plate(colonies,dims,direc,ff.chan,bIms,nIms, outfile);

    plate1.mm = 1;
    plate1.si = size(bIms{1});
    save([direc2 filesep outfile],'plate1','peaks','-append');  
disp('done');
