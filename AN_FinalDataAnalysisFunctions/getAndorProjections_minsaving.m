
% get max projections from minimum files Andor saving, merge all time
% groups

direc = 'E:\allSortingData\2018-04-18-liveSorting70to30CFPpluriS4diff\LiveSort_7030_cfpPluri_Smad4diff_20180418_125959 PM';

direc2save = 'E:\allSortingData\2018-04-18-liveSorting70to30CFPpluriS4diff\MIP_alltimes';
ff = readAndorDirectory(direc);
%ff1 = readAndorDirectory(direc2); 
for pos = ff.p(1)%:ff.p(end);
time = [0 1 2];
z = [];
for xx =1%;
filename1 = getAndorFileName(ff,pos,time(1),z,ff.w(xx));
filename2 = getAndorFileName(ff,pos,time(2),z,ff.w(xx));
filename3 = getAndorFileName(ff,pos,time(3),z,ff.w(xx));

r1 = bfGetReader(filename1);
r2 = bfGetReader(filename2);
%r3 = bfGetReader(filename3);


nz=r1.getSizeZ; % number of Z planes
tpts1 = r1.getSizeT;
tpts2 = r2.getSizeT;
tpts3 = 0;
%tpvect = 1:nz:(tpts1+tpts2);
multitp_nuc = [];
qq = 0;
for ii=1:(tpts1+tpts2+tpts3) % loop over total number of time points
    clear multitp_nuc
    reader = r1;
     if ii>tpts1
         reader = r2;
         qq = tpts1;
     end
     if ii>tpts1+tpts2
         reader = r3;
         qq = tpts1+tpts2;
     end
     multitp_nuc = bfMaxIntensity(reader,(ii-qq),1);
     %figure, imshow(multitp_nuc,[]);
     if (pos)<10
         imwrite(multitp_nuc,[direc2save '/' 'sortandregister_21to43hrs_MIP_f000' num2str(pos) '_w000' num2str(ff.w(xx)) '.tif'],'writemode','append');
     end
     if (pos)>=10
         imwrite(multitp_nuc,[direc2save '/' 'sortandregister_21to43hrs_MIP_f00' num2str(pos) '_w000' num2str(ff.w(xx)) '.tif'],'writemode','append');
     end
     disp(['saved projection for time point' num2str(ii) 'channel' num2str(ff.w(xx)) 'position' num2str(pos)]);
end
end
end


