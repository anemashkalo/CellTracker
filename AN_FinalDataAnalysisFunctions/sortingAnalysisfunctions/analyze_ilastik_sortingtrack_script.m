%% extract tracks from Ilastik file
% here fully use Ilastik Automatic tracking 
pos = 4;%1
chan = 1;% ff.w(chan): 0 - cfp cells; 1- nuc marker of other cell type
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
delta_t = 15;%15
datadir = 'E:\allSortingData\last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting';%last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting
ff_tmp = readAndorDirectory(datadir);
if ff_tmp.p(pos)<10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_last80tpts_MIP_f000' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];%SortingGFPS4cellspluri70to30_MIP_80tpts   SortingGFPS4cellspluri70to30_MIP_last80tpts
end
if ff_tmp.p(pos)>=10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_last80tpts_MIP_f00' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];
end
il_tracks
[coordintime] = process_ilastik_trackAN(il_tracks);
%% get the stats from the tracked movie and from the coresponding untracked cells in the same colony
ilastikprob = 1;
%paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
%global setUserParam 
direc_untracked ='E:\allSortingData\last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting';
ifile_untracked =...
'E:\allSortingData\last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting\SortingGFPS4cellspluri70to30_MIP_last80tpts_MIP_f0009_w0001_Probabilities.h5';%SortingGFPS4cellspluri70to30_MIP_80tpts_MIP_f0012_w0000_Probabilities.h5
chan_untracked = 1;% 0 - cfp cells; 1- nuc marker of other cell type
direc_tracked  =...
'E:\allSortingData\last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting';
chan_tracked = 0;
ifile_tracked =...
'E:\allSortingData\last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting\SortingGFPS4cellspluri70to30_MIP_last80tpts_MIP_f0009_w0000_Probabilities.h5';
[untrackedstats,trackedstats,img_untracked,img_tracked]=get_celltypes_statsimages(direc_untracked,ifile_untracked,direc_tracked,ifile_tracked,paramfile,pos,chan_tracked,chan_untracked,ilastikprob);
% TODO: save this data into file


%% analyze the motion of trackIDs selected using the GUI
% TODO: add the sarting ratio of the cell types
% the sorting happens faster in smaler colonies (f0000 vs f0001 ?)
close all
%untrackedstats,trackedstats,img_untracked,img_tracked
load('SelectedTrackIDsData_f9aftsort_chan_w0000.mat','goodtracks');%SelectedTrackIDsData_f12_chan_w0001
goodtracks(isnan(goodtracks)) = [];
goodtracks = unique(goodtracks);
local_neighbors = struct;
msd_all = struct;
frac_init = zeros(size(goodtracks,2),1);
frac_mean = zeros(size(goodtracks,2),1);
DiffC = zeros(size(goodtracks,2),1);
counter = 0;
toplot =1;%
toplot2 =0;
x1 = 8;
x2 = 55; % ideally this is the full track length, but if the tracks were observed to be ...
% good only untill certain time point, then that is the value of x2
allcells_speed_inT= struct;
allcells_local = struct;
total_disp_pertrack = zeros(size(goodtracks,2),1);
meanspeed_pertrack = zeros(size(goodtracks,2),1);
matlabtracking = 0;% was it ilastik-generated track or matlab-generated

for ii=1:size(goodtracks,2)
  %  if size(coordintime(ii).dat,1)>shortesttrack
    counter = counter+1;
    trackID =goodtracks(ii);    
        delta_t = 15;
        paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
        [samecell_neighbors,othercell_neighbors,same_neighbor_abs,other_neighbor_abs,curr_cell_speed,curr_cell_displ,global_neighborhood]=...
            getDynamicLocalNeighborsAN(trackID,paramfile,coordintime,matlabtracking,delta_t,untrackedstats,trackedstats,img_untracked,img_tracked,toplot,x2);%untrackedstats,trackedstats,img_untracked,img_tracked
        allcells_speed_inT(ii).dat = curr_cell_speed;
        allcells_disp_inT(ii).dat = curr_cell_displ;
        allcells_local(ii).same = samecell_neighbors;
        allcells_local(ii).other =othercell_neighbors;
        local_neighbors(ii).dat = othercell_neighbors;%samecell_neighbors
        [msd,slope_estimate,D,diff_coeff,mean_disp,total_disp,mean_speed,validTracks,MD,mspeed,TP,fit_out,totalD] = ...
        getCellmovement_params_IlastikTrack(coordintime,delta_t,trackID,x1,paramfile,othercell_neighbors,toplot2,x2);
        msd_all(ii).dat = msd(goodtracks(ii)).dat;
        total_disp_pertrack(ii,1) = totalD;
        meanspeed_pertrack(ii,1) = mspeed;
        if toplot2 ==1
            figure(1),hold on
            h = figure(1);
            h.Colormap = jet;
            caxis([0 1]);
            ylim([0 15000])
        end
        frac_init(ii,1) = (mean(local_neighbors(ii).dat(1:x1)));%othercell_neighbors(1:x1)1:x1
        frac_mean(ii,1) = (mean(local_neighbors(ii).dat(1:x2)));%isfinite
        DiffC(ii,1) = diff_coeff;
    %end
end
totalDispl = total_disp_pertrack;
[r,~]=find(isnan(frac_init));
frac_init(r) = [];
DiffC(r) = [];
total_disp_pertrack(r) = [];
cc = corrcoef(frac_init,DiffC);
cc = cc(1,2);
cc_disp = corrcoef(frac_init,total_disp_pertrack);%frac_mean
cc_disp = cc_disp(1,2);

figure(4), plot((total_disp_pertrack),frac_init,'p','Markersize',14,'MarkerFaceColor','r','MarkerEdgeColor','c');box on
h4 = figure(4);h4.Colormap = jet;ylim([0 1.1]);ylabel('Starting Fraction of other cell type in the neighborhood');xlabel('Total cell displacement per track, um');
title(['Initial fraction:first ' num2str(x1*delta_t/60) 'hrs of movie; Correlation coefficient ' num2str(cc_disp) ]);

figure(10), plot(frac_init,round(DiffC),'p','Markersize',13,'MarkerEdgeColor','k','MarkerFaceColor','g');box on;  hold on
xlabel(['Mean fraction of other cell type neighbors during ' num2str(x1*delta_t/60) 'hrs']);
title(['Total ' num2str(counter) ' uninterrupted tracks; Correlation coefficient ' num2str(cc) ]);
xlim([0 1.05]); 
ylim([0 round(max(DiffC)+5)]);
ylabel(['Estimated D [um^2/hr] from initial ~ ' num2str(x1*delta_t/60) '-hrs long motion']);

global_nbrhd = cat(1,global_neighborhood.frame);

figure(15), histogram(total_disp_pertrack,'normalization','probability','binwidth',15,'FaceColor','m');box on
ylim([0 1]);xlabel('Total cell displacement during tracked motion,um'); ylabel('Frequency')
%save(['qMotion_f000' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_time_' num2str(80*delta_t/60) 'hrs-' num2str((80+x2)*delta_t/60) 'hrs.mat' ],'x1','delta_t','DiffC','frac_init','allcells_speed_inT','counter','frac_mean','meanspeed_pertrack','x2','allcells_local','totalDispl','untrackedstats','trackedstats','global_nbrhd','allcells_disp_inT','same_neighbor_abs','other_neighbor_abs');

%save(['qMotion_f000' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_time_' num2str(1*delta_t/60) 'hrs-' num2str((x2)*delta_t/60) 'hrs.mat' ],'x1','delta_t','DiffC','frac_init','allcells_speed_inT','counter','frac_mean','meanspeed_pertrack','x2','allcells_local','totalDispl','untrackedstats','trackedstats','global_nbrhd','allcells_disp_inT','same_neighbor_abs','other_neighbor_abs');

[b,~]=find(isnan(frac_mean));
frac_mean(b) = [];
totalDispl(b) = [];
cc_disp2 = corrcoef(frac_mean,totalDispl);%frac_mean
cc_disp2 = cc_disp2(1,2);
figure(6), plot((totalDispl),frac_mean,'p','Markersize',14,'MarkerFaceColor','m','MarkerEdgeColor','b');box on
h6 = figure(6);h6.Colormap = jet;ylim([0 1.1]);ylabel('Fraction of other cell type in the neighborhood');xlabel('Total cell displacement per track, um');
title(['Fraction:mean during ' num2str(x2*delta_t/60) 'hrs of movie; Correlation coefficient ' num2str(cc_disp2) ]);

figure(7), histogram(meanspeed_pertrack,'normalization','probability','binwidth',3,'FaceColor','b');box on
ylim([0 1]);xlabel('Mean cell speed during tracked motion,um/hr'); ylabel('Frequency')

figure(23),plot((1:size(global_nbrhd,1))*delta_t/60,global_nbrhd,'b*'); box on
if max(global_nbrhd)<=1
title('Actual Ratio of cell types (CFP:PLURI)');ylim([0 1]);
else
title('Actual Ratio of cell types (PLURI:CFP)');ylim([0 2]);
end
xlabel('Time, hrs'); ylabel('Ratio');
%% plot diffusion coefficient distributions and mean speeds per track in time for two sorting cell types
close all
bnwdth =15;
pluriD = struct;
meanDisp_pluri = struct;
tracked_time_pluri= struct;
frac_initP = struct;
pluriV  = struct;
frac_meanP = struct;
colormap = jet;
matfile_str_w1 ={'qMotion_f0009_w0001_time_20hrs-33.75hrs.mat','qMotion_f00012_w0001_time_20hrs-32.5hrs.mat'};
colszR =100;
for q=1:size(matfile_str_w1,2)
load(matfile_str_w1{q},'x1','delta_t','DiffC','frac_init','allcells_speed_inT','counter','frac_mean','meanspeed_pertrack','x2','allcells_local','totalDispl','untrackedstats','trackedstats','global_nbrhd','allcells_disp_inT');%,'same_neighbor_abs','other_neighbor_abs'
pluriD(q).dat = DiffC;
meanDisp_pluri(q).dat = totalDispl;
tracked_time_pluri(q).t = x2;
frac_initP(q).dat = frac_init;
frac_meanP(q).dat = frac_mean;
pluriV(q).dat = meanspeed_pertrack;
end

pluriD_all = cat(1,pluriD.dat);
frac_initP_all = cat(1,frac_initP.dat);
meanDisp_pluri_all = cat(1,meanDisp_pluri.dat);
pluriV_all = cat(1,pluriV.dat);
frac_meanP_all=cat(1,frac_meanP.dat);
sz1=size(pluriD_all,1) ;

disp(['mean D pluri cells:' num2str(mean(pluriD_all)) 'std dev:' num2str(std(pluriD_all))]);
disp(['meandisplacement pluri cells:' num2str(mean(meanDisp_pluri_all)) 'std dev:' num2str(std(meanDisp_pluri_all)) ]);
disp(['mean speed pluri cells:' num2str(mean(pluriV_all)) 'std dev:' num2str(std(pluriV_all)) ]);
finDiffCoef =[colszR mean(pluriD_all) std(pluriD_all)];
finTDisplacement=[colszR mean(meanDisp_pluri_all)  std(meanDisp_pluri_all)];
finMeanVel=[colszR mean(pluriV_all)  std(pluriV_all)];

%save(['stats_pluri_colR' num2str(colszR) 'um.mat'],'finDiffCoef','finTDisplacement','finMeanVel','frac_initP_all');
figure(1),histogram(pluriD_all,'binwidth',bnwdth,'normalization','probability','FaceColor','r');hold on%,'binwidth',8
matfile_str_w0 ={'qMotion_f0009_w0000_time_20hrs-33.75hrs.mat','qMotion_f00012_w0000_time_20hrs-33.75hrs.mat'};
diffD = struct;
meanDisp_diff = struct;
tracked_time_diff= struct;
frac_initD = struct;
diffV  = struct;
frac_meanD = struct;
for q=1:size(matfile_str_w0,2)
load(matfile_str_w0{q},'x1','delta_t','DiffC','frac_init','allcells_speed_inT','counter','frac_mean','meanspeed_pertrack','x2','allcells_local','totalDispl','untrackedstats','trackedstats','global_nbrhd','allcells_disp_inT');%,'same_neighbor_abs','other_neighbor_abs'
diffD(q).dat = DiffC;
meanDisp_diff(q).dat = totalDispl;
tracked_time_diff(q).t = x2;
frac_initD(q).dat = frac_init;
frac_meanD(q).dat = frac_mean;
diffV(q).dat = meanspeed_pertrack;
figure(50),plot(((1:size(global_nbrhd,1))+80)*delta_t/60,global_nbrhd,'*','Color',colormap(randi(size(colormap,1)),:));hold on;box on
title('Actual Ratio of cell types (CFP:PLURI)');ylim([0 1]);
xlabel('Time, hrs');ylabel('Ratio, per frame,not colony');
leg_str{q} = matfile_str_w0{q}(9:14);
end
figure(50), legend(leg_str);
diffD_all = cat(1,diffD.dat);
sz2=size(diffD_all,1) ; 
frac_initD_all = cat(1,frac_initD.dat);
meanDisp_diff_all = cat(1,meanDisp_diff.dat);
diffV_all = cat(1,diffV.dat);
frac_meanD_all = cat(1,frac_meanD.dat);
disp(['meanD diff cells:' num2str(mean(diffD_all)) 'std dev:' num2str(std(diffD_all)) ]);
disp(['meandisplacement diff cells:' num2str(mean(meanDisp_diff_all)) 'std dev:' num2str(std(meanDisp_diff_all)) ]);
disp(['mean speed diff cells:' num2str(mean(diffV_all)) 'std dev:' num2str(std(diffV_all)) ]);
finDiffCoef =[colszR mean(diffD_all) std(diffD_all)];
finTDisplacement=[colszR mean(meanDisp_diff_all)  std(meanDisp_diff_all)];
finMeanVel=[colszR mean(diffV_all)  std(diffV_all)];

%save(['stats_cfp_colR' num2str(colszR) 'um.mat'],'finDiffCoef','finTDisplacement','finMeanVel','frac_initD_all');
figure(1),histogram(diffD_all,'binwidth',bnwdth,'normalization','probability','FaceColor','c');
ylim([0 1]);
ylabel('Frequency');
xlabel('Estimated D [um^2/hr]');
legend(['PLURI cells (' num2str(sz1) ') cells' ],['PRE-DIFF cells (' num2str(sz2) ') cells' ]);
title('Both inner and edge cells are considered')

figure(20), histogram(meanDisp_pluri_all,'FaceColor','r','normalization','probability','binwidth',10);hold on
histogram(meanDisp_diff_all,'FaceColor','c','normalization','probability','binwidth',10);ylim([0 0.6]);
title('Mean total displacement during tracked motion')
xlabel('Total displacement,um');ylabel('Frequency')
legend(['PLURI cells (' num2str(sz1) ') cells' ],['PRE-DIFF cells (' num2str(sz2) ') cells' ]);

% compare the D vs neighborhood dependence for cell types
figure(10), plot(frac_initP_all,pluriD_all,'p','Markersize',13,'MarkerEdgeColor','k','MarkerFaceColor','r');  hold on
plot(frac_initD_all,diffD_all,'p','Markersize',13,'MarkerEdgeColor','k','MarkerFaceColor','c');box on;
legend(['PLURI cells (' num2str(sz1) ') cells' ],['PRE-DIFF cells (' num2str(sz2) ') cells' ]);
xlabel(['Mean fraction of other cell type neighbors during ' num2str(x1*delta_t/60) 'hrs']);
xlim([0 1.05]); 
ylim([0 round(max(diffD_all)+5)]);
ylabel(['Estimated D [um^2/hr] from initial ~ ' num2str(x1*delta_t/60) '-hrs long motion']);
% [r,~]=find(isnan(frac_init));
% frac_init(r) = [];
% DiffC(r) = [];
cc = corrcoef((frac_initP_all),(pluriD_all));
cc1 = cc(1,2);
cc = corrcoef(frac_initD_all,(diffD_all));
cc2 = cc(1,2);
title(['Correlation coefficient for pluri cells ' num2str(cc1) ';  for CFP cells ' num2str(cc2) ]);
% compare distributions (no statistical toolbox for this)
%h = kstest2(round(pluriD),round(diffD));
%h = kstest2(x1,x2) returns a test decision for the null hypothesis that the data in vectors x1 and x2 are from the same continuous distribution, using the two-sample Kolmogorov-Smirnov test. The alternative hypothesis is that x1 and x2 are from different continuous distributions. The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
%velocities:
figure(2),histogram(pluriV_all,'FaceColor','r','binwidth',3,'normalization','probability');hold on
histogram(diffV_all,'FaceColor','c','binwidth',3,'normalization','probability');box on
ylabel('Frequency');ylim([0 1]);
legend(['PLURI cells (' num2str(sz1) ') cells' ],['PRE-DIFF cells (' num2str(sz2) ') cells' ]);
xlabel('Mean cell speed (averaged over track length), um/hr');
title('Average taken over track length');

figure(56),plot(meanDisp_diff_all,frac_meanD_all,'cp','Markersize',12,'LineWidth',1.2);hold on ;
plot(meanDisp_pluri_all,frac_meanP_all,'rp','Markersize',12,'LineWidth',1.2)
box on%diffV_all  meanDisp_diff_all  %pluriV_all  meanDisp_pluri_all
ylabel('Mean fraction of unlike cells throughout tracked motion');ylim([0 1]);
xlabel('Total cell displacement, um');
legend(['PLURI cells (' num2str(sz1) ') cells' ],['PRE-DIFF cells (' num2str(sz2) ') cells' ]);
%%
%  todo2: absolute number of neighboring cells as a measure of ecad
%  molecules levels vs cell speed/displ; not saved correctly above



%% todo1: change to displacement vs local neighborhood
close all
load('qMotion_f0009_w0001_time_20hrs-31.5hrs.mat','DiffC','frac_init','allcells_speed_inT','counter','x2','allcells_local');ylim([0 0.8]);
pluriV = allcells_speed_inT;
pluri_local = allcells_local;
for jj=1:size(pluriV,2)
time = 1:x2;%size(pluriV(jj).dat,1);
figure(1),scatter(pluriV(jj).dat(1:x2),pluri_local(jj).other(1:x2),[],time);hold on
end
szp = x2;
h = figure(1);
h.Colormap = jet;
caxis([0 max(szp)]);
ylim([0 1]);box on
ylabel('Fraction of other type neighbors');
xlabel('Speed, um/hr');colorbar
%% comparison of parameters of motion as var colony sizes
close all
D_pluri = [51.0 28.7 27.2];
colszP = [100 400 650 ];
D_diff= [38.4 36.1 32.8];
colszD = [100 400 650 ];
figure(1),plot(colszP,D_pluri,'rp','Markersize',13,'LineWidth',1.5);hold on
plot(colszD,D_diff,'bp','Markersize',18,'LineWidth',1.5);box on;hold on

D_pluri = [20.0];
colszP = [100];
D_diff= [22.6];
colszD = [100 ];
plot(colszP,D_pluri,'r*','Markersize',13,'LineWidth',1.5);hold on
plot(colszD,D_diff,'b*','Markersize',18,'LineWidth',1.5);box on;hold on
title('Mean D coeff,(each data point is mean over ~ 20-40 cells)');
xlabel('Colony size, R, um');
ylabel('D estimate, um^2/hr');
xlim([10 1000]);
ylim([15 65]);
legend('Pluri cells','diff cells')

% redo here, get the acutal mean values
cc_DvsLocal_pluri=[0.04 0.45 -0.07];
DvsLocal_sz1 = [100 400 650 ];
cc_DvsLocal_diff=[0.38 0.54 0.2];
DvsLocal_sz2 = [100 400 650];
figure(2),plot(DvsLocal_sz1,cc_DvsLocal_pluri,'rp','Markersize',13,'LineWidth',1.5);hold on
plot(DvsLocal_sz2,cc_DvsLocal_diff,'bp','Markersize',18,'LineWidth',0.7);box on
title('Corr. coeff. btw D and fraction of unlike cells during initial motion');
legend('Pluri cells','diff cells')
xlabel('Colony size, R, um');
ylabel('Correlation coefficient');
xlim([10 1000]);
ylim([-0.2 0.8]);



mean_cell_speed_diff = [13.4 13.7 13.4];
v_sz1 = [100 400 650];
mean_cell_speed_pluri = [16.9 13.8 14.9];
v_sz2 = [100 400 650];
figure(3),plot(v_sz1,mean_cell_speed_diff,'bp','Markersize',20,'LineWidth',1.2);hold on
plot(v_sz2,mean_cell_speed_pluri,'rp','Markersize',15,'LineWidth',1.2);box on;hold on
mean_cell_speed_diff = [9.8];
v_sz1 = [100 ];
mean_cell_speed_pluri = [12.4];
v_sz2 = [100];
plot(v_sz1,mean_cell_speed_diff,'b*','Markersize',20,'LineWidth',1.2);hold on
plot(v_sz2,mean_cell_speed_pluri,'r*','Markersize',15,'LineWidth',1.2);box on;hold on
title('Mean cell speed during 11 hours,(each data point is mean over ~ 20-40 cells)');
xlabel('Colony size, R, um');
ylabel('Speed, um/hr');
xlim([10 1000]);
ylim([8 19]);
legend('diff cells','Pluri cells');

total_disp_diff = [ 31.3 45.6 33.1];
v_sz1 = [100 400 650];
total_disp_pluri = [40.8 24.2 42.6];
v_sz2 = [100 400 650];
figure(4),plot(v_sz2,total_disp_pluri,'rp','Markersize',11,'LineWidth',1.2);hold on
plot(v_sz1,total_disp_diff,'bp','Markersize',20,'LineWidth',1.2);box on ;hold on
total_disp_diff = [34.2];
v_sz1 = [100 ];
total_disp_pluri = [32.9];
v_sz2 = [100];
figure(4),plot(v_sz2,total_disp_pluri,'r*','Markersize',11,'LineWidth',1.2);hold on
plot(v_sz1,total_disp_diff,'b*','Markersize',20,'LineWidth',1.2);box on ;hold on
title('Mean total displacement during 11 hrs (each point is mean over ~ 20-40 cells)');
legend('Pluri cells','diff cells')
xlabel('Colony size, R, um');
ylabel('Cell displacement, um');
xlim([10 1000]);
ylim([20 55]);



