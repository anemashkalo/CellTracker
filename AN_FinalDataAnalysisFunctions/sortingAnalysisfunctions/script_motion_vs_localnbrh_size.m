%1. whrte a function where the local neighborhood size is varied
%2. motion params in several dts
%3. D coeff in windows of several hours
%4. Smad4 vs nighborhood (colorcoded intime)
%5. error bars everywhere
%6. spatial correlation btw displacements of cells of diff type
%% get the tracking data from ilastik file 
pos = 2;%1
chan = 1;% ff.w(chan): 0 - cfp cells; 1- nuc marker of other cell type
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
delta_t = 15;%15
datadir = 'E:\allSortingData\80_tp_fortracking_2017-07-14-S4cfpDataset_sorting';%last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting
ff_tmp = readAndorDirectory(datadir);
if ff_tmp.p(pos)<10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_80tpts_MIP_f000' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];%SortingGFPS4cellspluri70to30_MIP_80tpts   SortingGFPS4cellspluri70to30_MIP_last80tpts
end
if ff_tmp.p(pos)>=10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_80tpts_MIP_last80tpts_MIP_f00' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];
end
il_tracks
[coordintime] = process_ilastik_trackAN(il_tracks);
%% get images
pos=2;
chan_untracked=1;
direc_untracked='E:\allSortingData\80_tp_fortracking_2017-07-14-S4cfpDataset_sorting';
ff1 = readAndorDirectory(direc_untracked);
nucmoviefile1 = getAndorFileName(ff1,ff1.p(pos),[],[],chan_untracked);%getAndorFileName(files,pos,time,z,w)
img_untracked = bfopen(nucmoviefile1);

direc_tracked='E:\allSortingData\80_tp_fortracking_2017-07-14-S4cfpDataset_sorting';
chan_tracked=0;
ff2 = readAndorDirectory(direc_tracked);
nucmoviefile = getAndorFileName(ff2,ff2.p(pos),[],[],chan_tracked);%getAndorFileName(files,pos,time,z,w)
% nreader = bfGetReader(nucmoviefile);
% nt = nreader.getSizeT;
img_tracked = bfopen(nucmoviefile);
%%
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
matfile_1 = 'qMotion_f0001_w0000_time_0.25hrs-11.25hrs.mat';
matfile_2 = 'SelectedTrackIDsData_f1_chan_w0000.mat';
load(matfile_1);
load(matfile_2);
var_sz = [35 45 55 65 85 95];% 
run(paramfile)
allcells_disp_inT=struct;
heteroneighbors_inT=struct;
homoneighbors_inT = struct;
global userParam
N= 8;
str1 = 'CFP cells';
mean_ccother=[];
err_ccother=[];
mean_ccsame=[];
err_ccsame =[];
cc2=[];
cc1=[];
for jj=1:size(var_sz,2)
    
for w = 1:size(goodtracks,2)
    trackID=goodtracks(w);
[samecell_neighbors,othercell_neighbors,curr_cell_speed,curr_cell_displ]= ...
    vary_localnbr_range(trackID,paramfile,coordintime,delta_t,untrackedstats,trackedstats,img_untracked,img_tracked,x2,var_sz(jj));
allcells_disp_inT(w).dat = curr_cell_displ;
heteroneighbors_inT(w).dat = othercell_neighbors;
homoneighbors_inT(w).dat = samecell_neighbors;
d = allcells_disp_inT(w).dat;%allcells_disp_inT
nbr_other = heteroneighbors_inT(w).dat;
nbr_same= homoneighbors_inT(w).dat;
tmp1=[];
tmp2 = [];
tmp3 = [];
q = 1;
for k=1:N:(x2-N+1)
    tmp1(q) = mean(d(k:k+N-1));
    tmp2(q) = mean(nbr_other(k:k+N-1));
    tmp3(q) = mean(nbr_same(k:k+N-1));
    %disp([k k+N-1]);   
    q = q+1;
end
    cc_other = corrcoef(tmp1,tmp2);
    cc_same = corrcoef(tmp1,tmp3);
   if isfinite(cc_other(2)) && (isfinite(cc_same(2)) )
    cc1(w,1) = cc_other(2);
    cc2(w,1) = cc_same(2);
   end
end
mean_ccother(jj) = mean(cc1(cc1>0));%(cc2>0)
err_ccother(jj) = std(cc1(cc1>0));%
mean_ccsame(jj) = mean(cc2(cc2>0));%
err_ccsame(jj) = std(cc2(cc2>0));%

end

figure(2), errorbar(var_sz(isfinite(mean_ccother))*userParam.pxtomicron,mean_ccother(isfinite(mean_ccother)),err_ccother(isfinite(mean_ccother)),'p','Markersize',5,'MarkerFaceColor','b');hold on;
box on
xlabel('Local neighborhood size, um')
ylabel('Corr.coeff. btw cell displacement and unlike cell neighbors fraction');
title(['Data averaged every ' num2str(N*delta_t/60) 'hrs']);
ylim([0 1])
figure(3), errorbar(var_sz(isfinite(mean_ccsame))*userParam.pxtomicron,mean_ccsame(isfinite(mean_ccsame)),err_ccsame(isfinite(mean_ccsame)),'p','Markersize',5,'MarkerFaceColor','b');hold on;
box on
xlabel('Local neighborhood size, um')
ylabel('Corr.coeff. btw cell displacement and same-cell neighbors fraction');
title(['Data averaged every ' num2str(N*delta_t/60) 'hrs']);
ylim([0 1])

