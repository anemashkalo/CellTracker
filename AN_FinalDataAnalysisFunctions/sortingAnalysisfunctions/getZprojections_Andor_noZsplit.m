% direc = directory with raw images
direc = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-31-EcadReporterYFPwithbetacatRFP_liveimaging';
% direc2save = where to save the max. projections
direc2save = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-31-EcadReporterYFPwithbetacatRFP_liveimaging\MaxProjections';%\MaxProjections
ff = readAndorDirectory(direc);
delta_t = 30;% time interval between the frames in a given movie
totest = [ff.p];% which positions in the folder will be processed to make max projections 
tgroup = 0;
chan1 =[1];

for ii = 1:size(totest,2)
    pos = totest(ii);
%tpname = getAndorFileName(ff,pos,tgroup,[],ff.w(chan1));%chan
tpname = getAndorFileName(ff,pos,tgroup,[],[]);%chan
reader = bfGetReader(tpname);
nz=reader.getSizeZ;
nT = reader.getSizeT;
%nT = 80;
tpt1 = 1;
tpt2 = nT;
multitp_nuc = [];
for time=1:nT 
    for chan = 1:size(chan1,2)
        for jj = 1:nz
            iPlane=reader.getIndex(jj - 1, chan -1, time - 1) + 1;
            img_now=bfGetPlane(reader,iPlane);
            %figure(ii),imshow(img_now,[]);
            if jj == 1
                max_img =  img_now;
            else
                max_img = max(max_img,img_now);
            end
            multitp_nuc = (max_img);
        end
        
        if (pos)<10
            %ff.w(chan1)=chan-1;
            imwrite(multitp_nuc,[direc2save '\' ff.prefix '_hr' num2str(tpt1*delta_t/60) 'to' num2str(tpt2*delta_t/60) '_MIP_f000' num2str(pos) '_w000' num2str(chan-1) '.tif'],'writemode','append');
        end
        if (pos)>=10
            %ff.w(chan1)=chan-1;
            imwrite(multitp_nuc,[direc2save '\' ff.prefix '_hr' num2str(tpt1*delta_t/60) 'to' num2str(tpt2*delta_t/60) '_MIP_f00' num2str(pos) '_w000' num2str(chan-1) '.tif'],'writemode','append');
        end
        disp(['saved projection for time point' num2str(time) ' channel' num2str(chan-1) '  position' num2str(pos)]);
        
    end
end
end
 