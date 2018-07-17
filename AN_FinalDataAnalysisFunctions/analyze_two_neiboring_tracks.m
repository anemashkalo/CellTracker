%% get the tracking data
pos =1%1
chan = 1;% ff.w(chan): 0 - cfp cells; 1- nuc marker of other cell type
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
datadir = 'E:\allSortingData\80_tp_fortracking_2017-07-14-S4cfpDataset_sorting';%80_tp_fortracking_2017-07-14-S4cfpDataset_sorting last80_tp_geneexpression_2017-07-14-S4cfpDataset_sorting
cd(datadir)
ff_tmp = readAndorDirectory(datadir);
if ff_tmp.p(pos)<10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_80tpts_MIP_f000' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];%SortingGFPS4cellspluri70to30_MIP_80tpts   SortingGFPS4cellspluri70to30_MIP_last80tpts
end
if ff_tmp.p(pos)>=10
il_tracks = ['SortingGFPS4cellspluri70to30_MIP_80tpts_MIP_f00' num2str(ff_tmp.p(pos)) '_w000' num2str(ff_tmp.w(chan)) '_Tracking-Result.h5'];
end
disp(il_tracks);
[coordintime] = process_ilastik_trackAN(il_tracks);
%% script, get the stats of the two neighboring cells
close all 
clear all
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
motiondata_tracked = 'qMotion_f0000_w0000_time_0.25hrs-11.25hrs.mat';
load(motiondata_tracked,'trackedstats','untrackedstats','x2','delta_t','coordintime');
good_tracksfile = 'SelectedTrackIDsData_f0_chan_w0000.mat';
load(good_tracksfile);
trackID = goodtracks(9);%done:1,3,4,5,6,8,
toplot =1;
toshow = 1;
[sametype,othertype,tracked_coord,Displ_trackedcell]=get_local_cell_stats(coordintime,trackID,untrackedstats,trackedstats,x2,paramfile,delta_t,toplot,toshow);

%% test the output
close all;tp = 1;goodtracks;
figure(2), imshowpair(untrackedstats(tp).img,trackedstats(tp).img);hold on
    plot(tracked_coord(tp,1),tracked_coord(tp,2),'pc','Markersize',8,'MarkerFaceColor','c','MarkerEdgeColor','r');hold on
    if ~isempty(sametype(tp).dat)
    plot(sametype(tp).dat(:,1),sametype(tp).dat(:,2),'co','Markersize',8,'MarkerFaceColor','c');hold on
    %text(sametype(tp).dat(:,1),sametype(tp).dat(:,2),num2str(sametype(tp).dat(:,3)),'Color','k');hold on
    end
    if ~isempty(othertype(tp).dat)
    plot(othertype(tp).dat(:,1),othertype(tp).dat(:,2),'bo','Markersize',8,'MarkerFaceColor','b');hold on%,'LineWidth',0.8
   end
        title(['Time point  ' num2str(delta_t*tp/60) ' hrs; frame: ' num2str(tp) ]);
disp(['At first time point there are ' num2str(size(othertype(1).dat,1)) ' unlike type and '   num2str(size(sametype(1).dat,1)) '  like type neighbors' ])
title(['Time point  ' num2str(delta_t*tp/60) ' hrs; frame: ' num2str(tp) ]);
disp(['At first time point there are ' num2str(size(othertype(1).dat,1)) ' unlike type and '   num2str(size(sametype(1).dat,1)) '  like type neighbors' ])        
%% now analyze the two tracks (homotypic or heterotypic)

close all;run(paramfile);global userParam
xy_towatch = zeros(x2,3);
alldata = struct;
cellN =10;

samecelltype =1;
if samecelltype == 1
    alldata  = sametype;
else
    alldata  = othertype;
end
xy_towatch(1,1:2)= alldata(1).dat(cellN,1:2);% 
xy_towatch(1,3)= 1;
closest_cell=[];
for jj=1:size(alldata,2)-1 
    if ~isempty(alldata(jj+1).dat)
   tmp = ipdm(xy_towatch(jj,1:2),alldata(jj+1).dat(:,1:2),'Result','Structure','Subset','NearestNeighbor');%
   if tmp.distance <=  userParam.maxdist_tomove % make sure that a different cell was not picked up,if so,don't store or plot
   closest_cell = alldata(jj+1).dat(tmp.columnindex,1:2);
   % add the pixel overlap tooalldata(jj+1).dat(:,1:2)
   xy_towatch(jj+1,1:2) = closest_cell;
   xy_towatch(jj+1,3) = jj+1;    
   figure(3), imshowpair(untrackedstats(jj).img,trackedstats(jj).img);hold on  
plot(tracked_coord(jj,1),tracked_coord(jj,2),'pc','Markersize',8,'MarkerFaceColor','c','MarkerEdgeColor','r');hold on
if samecelltype == 1
plot(xy_towatch(jj,1),xy_towatch(jj,2),'bo','Markersize',8,'MarkerFaceColor','c');hold on
else
plot(xy_towatch(jj,1),xy_towatch(jj,2),'bo','Markersize',7,'MarkerFaceColor','r');hold on
end
title(['Time point: ' num2str(delta_t*jj/60) ' hrs;frame(' num2str(jj) ')Cells considered within ~ ' num2str(round(userParam.pxtomicron*userParam.local_sz)) 'um of tracked cell' ]);
   end
    end   
end

%%
[d_btw]=get_dist_btw_cellcenters(xy_towatch,tracked_coord,paramfile,delta_t);
figure(5), plot(d_btw(:,2),d_btw(:,1),'*m');hold on;box on
%errorbar(d_btw(1:4:end,2),mean(d_btw(1:4:end,1))*ones(size(d_btw(1:4:end,2),1),1),std(d_btw(1:4:end,1))*ones(size(d_btw(1:4:end,2),1),1),'pb');
%errorbar(round(d_btw(end,2)/2),mean(d_btw(:,1)),std(d_btw(:,1)),'pb');
%d_btw_allunlike = struct;
h5 = figure(5);h5.Colormap = prism;
h5.CurrentAxes.FontSize = 12;
h5.CurrentAxes.LineWidth = 2;
xlabel('Time, hrs')
ylabel('Distance between centroids of two cells, um')
if samecelltype == 1
title(['Tracked cell and nearest like cell, mean D = ' num2str(mean(d_btw(:,1))) 'um'])
else
title(['Tracked cell and nearest unlike cell, mean D = ' num2str(mean(d_btw(:,1))) 'um' ])
end
ylim([0 max(d_btw(:,1))+2])
%% save like cells data
load('distance_btw_Like_cells_f0w0');
  ii = 23;
 d_btw_all_like(ii).dat =d_btw(:,:);%1:40
 if ii == 1
 save('distance_btw_Like_cells_f0w0','d_btw_all_like','samecelltype');%
 else
      save('distance_btw_Like_cells_f0w0','d_btw_all_like','samecelltype','-append');%

 end
%% save unlike cells
load('distance_btw_unlike_cells_f0w0');
  ii =14;
  d_btw_allunlike(ii).dat = d_btw(:,:);%1:18
  if ii == 1
  save('distance_btw_unlike_cells_f0w0','d_btw_allunlike','samecelltype')%
  else
        save('distance_btw_unlike_cells_f0w0','d_btw_allunlike','samecelltype','-append')%

  end
%% plot unlike cells data
close all
ff = dir('.');
alldat = [];
for ii =1:size(ff)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'distance_btw_unlike'))
        disp(['loading: ' num2str(ff(ii).name)])
        load(ff(ii).name);
        %disp(size(d_btw_allunlike))
        alldat = [d_btw_allunlike alldat];
    end
end
c = 0;
init_separation_d = 30;
meandat = [];
for jj=1:size(alldat,2)
 if ~isempty(alldat(jj).dat) && (alldat(jj).dat(1,1)<=init_separation_d )%size(alldat(jj).dat(:,2),1)
c = c+1;
     figure(5), plot(alldat(jj).dat(:,2),alldat(jj).dat(:,1),'-*');hold on;box on
   for kk=1:size(alldat(jj).dat,1)
       meandat(kk,jj) = alldat(jj).dat(kk,1);
   end
 end
    end
ylabel('Distance between centroids of two cells, um')
title(['Unlike type cells start less than ' num2str(init_separation_d) 'um of each other'])
xlabel('time,hrs')
delta_t = 15;
for i=1:size(meandat,1)
figure(5), hold on
plot(i*(delta_t)/60,mean(nonzeros(meandat(i,:))),'kp','MarkerSize',13,'MarkerFaceColor','r');hold on
ylim([5 70])

figure(6),errorbar(i*(delta_t)/60,mean(nonzeros(meandat(i,:))),std(nonzeros(meandat(i,:)))./power(size(nonzeros(meandat(i,:)),1),0.5),'p','MarkerSize',13,'MarkerFaceColor','r','MarkerEdgeColor','k');box on;hold on

end
text(1,41,['total cells: ' num2str(c)],'Fontsize',12);

ylabel('Distance btw centroids of 2 cells (um),error bar = SEM')
title(['Unlike type cells start less than ' num2str(init_separation_d) 'um of each other'])
xlabel('time,hrs')
ylim([15 42])

%% plot like cells data
close all
ff = dir('.');
alldat = [];
for ii =1:size(ff)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'distance_btw_Like'))
        disp(['loading: ' num2str(ff(ii).name)])
        load(ff(ii).name);
        %disp(size(d_btw_all_like))
        alldat = [d_btw_all_like alldat];
    end
end
c = 0;
init_separation_d = 30;
meandat = [];

for jj=1:size(alldat,2)
 if ~isempty(alldat(jj).dat) && (alldat(jj).dat(1,1)>=init_separation_d)%size(alldat(jj).dat(:,2),1)
 c = c+1;
     figure(5), plot(alldat(jj).dat(1:end,2),alldat(jj).dat(1:end,1),'-*');hold on;box on
     for kk=1:size(alldat(jj).dat,1)
         meandat(kk,jj) = alldat(jj).dat(kk,1);
     end
 end
end
ylabel('Distance between centroids of two cells, um')
title(['Same type cells start more than ' num2str(init_separation_d) 'um of each other'])
xlabel('time,hrs')
delta_t = 15;
for i=1:size(meandat,1)
figure(5), hold on
plot(i*(delta_t)/60,mean(nonzeros(meandat(i,:))),'p','MarkerSize',13,'MarkerFaceColor','r','MarkerEdgeColor','k');hold on
ylim([15 70])

figure(6),errorbar(i*(delta_t)/60,mean(nonzeros(meandat(i,:))),std(nonzeros(meandat(i,:)))./power(size(nonzeros(meandat(i,:)),1),0.5),'p','MarkerSize',13,'MarkerFaceColor','c','MarkerEdgeColor','k');box on;hold on%
end
text(1,44,['total cells: ' num2str(c)],'Fontsize',12);

ylabel('Distance btw centroids of 2 cells (um),error bar = SEM')
title(['Same type cells start more than ' num2str(init_separation_d) 'um of each other'])
xlabel('time,hrs')
ylim([20 45])
%% calculate the displacement of each cell type and find the correlation
trackedcellD = cat(1,Displ_trackedcell.dat);%
[nearcellD,actual_times,dx_nearcell,dy_nearcell,dx_trackcell,dy_trackcell]=get_cell_displ(xy_towatch,paramfile,delta_t,tracked_coord);
figure(4), scatter(nearcellD,trackedcellD(1:actual_times-1),[],(1:actual_times-1),'LineWidth',2);box on
figure(7), scatter(cat(1,dx_nearcell.dat)./(delta_t/60),cat(1,dx_trackcell.dat)./(delta_t/60),[],(1:actual_times-1),'LineWidth',2);box on
xlabel('dx/dt of Immediate neighbor cell, um');
ylabel('dx/dt of Tracked cell, um');

h = figure(4);
h.Colormap = jet;
h.CurrentAxes.LineWidth  =3;
h.CurrentAxes.FontSize  =14;
xlim([0 max(cat(1,nearcellD,trackedcellD(1:actual_times-1)))]);
ylim([0 max(cat(1,nearcellD,trackedcellD(1:actual_times-1)))])
colorbar
caxis([1 actual_times-1])
xlabel('Immediate neighbor cell displacement, um')
ylabel('Tracked cell displacement, um')
cc =corrcoef(cat(1,dx_nearcell.dat),cat(1,dx_trackcell.dat))
if  samecelltype == 1
title(['Same cell type neighbor within' num2str(round(userParam.pxtomicron*userParam.local_sz)) 'um of tracked cell']);
else
  title(['Other cell type neighbor within ' num2str(round(userParam.pxtomicron*userParam.local_sz)) 'um of tracked cell']);
end
%% save mutual cell displacements
%mutual_displacement  =struct;
if samecelltype == 0
ii = 14;    
mutual_displacement(ii).tracked =trackedcellD(1:actual_times-1);
mutual_displacement(ii).neighbor =nearcellD;
%ii = ii+1;
%save('displacement_btw_unlike_cells','mutual_displacement','dx_nearcell','dy_nearcell','dx_trackcell','dy_trackcell','-append');%
else
ii = 15;    
mutual_displacement(ii).tracked =trackedcellD(1:actual_times-1);
mutual_displacement(ii).neighbor =nearcellD;
%save('displacement_btw_Like_cells','mutual_displacement','dx_nearcell','dy_nearcell','dx_trackcell','dy_trackcell','-append');% 
end
%% look an mutual displacement between starting cell pairs as like-like or unlike-like celll pair
counter = 0;
load('displacement_btw_Like_cells','mutual_displacement');%displacement_btw_Like_cells
for k=1:size(mutual_displacement,2)
    if ~isempty(mutual_displacement(k).tracked)
        counter = counter+1;
figure(9), scatter(mutual_displacement(k).tracked,mutual_displacement(k).neighbor,[],(1:size(mutual_displacement(k).tracked,1)));hold on
    end
    end