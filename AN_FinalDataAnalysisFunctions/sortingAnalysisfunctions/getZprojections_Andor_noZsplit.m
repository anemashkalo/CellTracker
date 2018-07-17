direc = 'E:\allSortingData\2017-06-29-livesorting_SDConfocalbetaCatcellspluri\livesortingBetacatpluri_20170704_20933 PM';
% dt = x mins 0.617284 um/pxl (20X)
direc2save = 'E:\allSortingData\Initial_tp_betaCatEsiCFP_separateimages';
ff = readAndorDirectory(direc);
delta_t = 12;
totest = [5];%%0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
tgroup = 0;
chan1 =3;
for ii = 1:size(totest,2)
    pos = totest(ii);
tpname = getAndorFileName(ff,pos,tgroup,[],ff.w(chan1));%chan
reader = bfGetReader(tpname);
nz=reader.getSizeZ;
nT = reader.getSizeT;
%nT = 80;
tpt1 = 50;
tpt2 = 50;
multitp_nuc = [];
for time=50%:nT  %1:nT%80:178
    chan = 1;
    for ii = 1:nz
        %time = 3;
        iPlane=reader.getIndex(ii - 1, chan -1, time - 1) + 1;
        img_now=bfGetPlane(reader,iPlane);
        %figure(ii),imshow(img_now,[]);
        if ii == 1
            max_img =  img_now;
        else
            max_img = max(max_img,img_now);
        end
        multitp_nuc = max_img;
    end
    if (pos)<10
         imwrite(multitp_nuc,[direc2save '\' ff.prefix '_hr' num2str(tpt1*delta_t/60) 'to' num2str(tpt2*delta_t/60) '_MIP_f000' num2str(pos) '_w000' num2str(ff.w(chan1)) '.tif'],'writemode','append');
    end
    if (pos)>=10
         imwrite(multitp_nuc,[direc2save '\' ff.prefix '_hr' num2str(tpt1*delta_t/60) 'to' num2str(tpt2*delta_t/60) '_MIP_f00' num2str(pos) '_w000' num2str(ff.w(chan1)) '.tif'],'writemode','append');
    end
          disp(['saved projection for time point' num2str(time) ' channel' num2str(ff.w(chan1)) '  position' num2str(pos)]);

end
end
 