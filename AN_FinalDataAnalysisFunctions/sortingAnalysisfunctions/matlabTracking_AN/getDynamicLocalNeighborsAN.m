% get the centroid coordinatesof the other cell type (non-CFP)
function [samecell_neighbors,othercell_neighbors,same_neighbor_abs,other_neighbor_abs,curr_cell_speed,curr_cell_displ,global_neighborhood]= getDynamicLocalNeighborsAN(trackID,paramfile,coordintime,matlabtracking,delta_t,pluristats,diffstats,img_pluri,img_diff,toplot,x2)
run(paramfile)
global userParam

for totrack =trackID 
tpts = size(coordintime(totrack).dat,1);%size(cfpstats,2)
time = 3;
if matlabtracking == 1
    time = 5;
    tpts = size(coordintime(totrack).dat,1);
end
colormap = prism;
v_2 = struct; % totla velocity at each time point
v = struct; %
othertype_neighbor = struct;
sametype_neighbor = struct;
fraction_same = struct;
fraction_other = struct;
global_neighborhood = struct;% what is the actual ratio of cell types (in the frame,not colony)
direction = struct;
displ = struct;
for jj = 1: (tpts-1)% time points    (size(cfpstats,2)-1)
tmp1 = round(cat(1,pluristats(jj).stats.Centroid));% centrois of all pluri cells at tp jj
area1 = cat(1,pluristats(jj).stats.Area);
allcells_type1 = cat(2,tmp1,area1);
tmp2 = round(cat(1,diffstats(jj).stats.Centroid));% centrois of all pluri cells at tp jj;
area2 = cat(1,diffstats(jj).stats.Area);
allcells_type2 = cat(2,tmp2,area2);
curr_track = coordintime(totrack).dat(jj,:);% coord of cell of other type, for which the neighborhood is quantified
rawimg1 = (img_pluri{1}{jj});% untracked
rawimg1 = imadjust(rawimg1,stretchlim(rawimg1));
rawimg = (img_diff{1}{jj});% tracked
rawimg = imadjust(rawimg,stretchlim(rawimg));
global_neighborhood(jj).frame = size(allcells_type2,1)/size(allcells_type1,1); % tracked to untracked
total_img = max(rawimg1,rawimg);
%at each time point find how many cells of each type are surrounding
% the given cell ( within the several cell radius)
allcells = cat(1,allcells_type1,allcells_type2); % all cells, including the coord of a current cell 
local_neighbors = ipdm(coordintime(totrack).dat(jj,1:2),allcells(:,1:2),'Result','Structure','Subset','Maximum','Limit',userParam.local_sz);% get all the cells, closer than local_sz
% v = power((vx^2+vy^2),0.5); vx = dx/dt;
tp1 =jj;
tp2 = jj+1;
% disp([tp1 tp2]);
dt = (tp2-tp1)*delta_t/60;
dx = (coordintime(totrack).dat(tp2,1)-coordintime(totrack).dat(tp1,1))*userParam.pxtomicron;
dy = (coordintime(totrack).dat(tp2,2)-coordintime(totrack).dat(tp1,2))*userParam.pxtomicron;
vx_2 = power(dx/dt,2);
vy_2 = power(dy/dt,2);
v_2(jj).total = power((vx_2+vy_2),0.5); 
v(jj).vx = dx/dt;
v(jj).vy = dy/dt;
displ(jj).dat = power(dx*dx+dy*dy,0.5);
% TODO: calculate the angle between the velocity vector of the tracked cell
% and the velocity vector of the cell of the opposite type (one from the
% neighborhood). Or do the angle btw the v(cfp cell) and the refernce direction +Ox, for
% example
if (v(jj).vy>=0) && (v(jj).vx>=0)
direction(jj).rad = atan(abs(v(jj).vy)/abs(v(jj).vx));% returns angle in radians  
direction(jj).theta = (direction(jj).rad)*180/3.14; % in degrees
end
%----------
if (v(jj).vy>=0) && (v(jj).vx<=0)
direction(jj).rad = 3.14-atan(abs(v(jj).vy)/abs(v(jj).vx));% returns angle in radians      
direction(jj).theta = 180-((direction(jj).rad)*180/3.14); % in degrees
end
if (v(jj).vy<=0) && (v(jj).vx>=0)
direction(jj).rad = (3.14*270)/180+atan(abs(v(jj).vy)/abs(v(jj).vx));% returns angle in radians
direction(jj).theta = 270+((direction(jj).rad)*180/3.14); % in degrees
end
if (v(jj).vy<=0) && (v(jj).vx<=0)
direction(jj).rad = 3.14+atan(abs(v(jj).vy)/abs(v(jj).vx));% returns angle in radians   
direction(jj).theta = 180+((direction(jj).rad)*180/3.14); % in degrees
end
%----------
% need to remove the closest cell from the list (since it's the tracked
% cell itself):
[~,c]=find(local_neighbors.distance==min(local_neighbors.distance));
local_neighbors.columnindex(c)=[];
% local_neighbors.columnindex - cells within the rad of local_sz
% now find which of these neighbors are cfp and wich are pluri
% then see if there is a correlaion between celltype1 motion (velocity of tracked cell,
% etc) and the number of neighbors of specific type
totest1 = allcells(local_neighbors.columnindex',:);
totest = totest1(totest1(:,3)>=userParam.arealow_localnbr,:);% filter for area  in the identified local neighborhood
nearest = size(totest,1);% all the cells within the local_sz pixels of the tracked cell
counter = 0;
for h=1:size(totest,1)
    % if one of the found local neighborhood cells are in the set of
    % celltype2, the nearest neighbor to it will be at zero distance (the cell is closest to itself);
    % counting how many of those, gives the number of neighbors of cell
    % type2, the rest (nearest-counter) is the other cell type
    tmp = ipdm(totest(h,1:2),allcells_type2(:,1:2),'Result','Structure','Subset','NearestNeighbor');%    
    if tmp.distance == 0
        counter = counter+1;        
    end
end

if (jj == 11) && (toplot == 1)%last tp: (tpts-1) or x2
%figure(jj),imshow(total_img,[500 1500]);hold on % 
figure(jj),imshowpair(rawimg1,rawimg);hold on % show both cell types onn one image
figure(jj),plot(coordintime(totrack).dat(jj,1),coordintime(totrack).dat(jj,2),'kp','MarkerFaceColor','y','MarkerSize',12,'LineWidth',1);hold on%colormap(randcolor,:)
figure(jj),plot(totest(:,1),totest(:,2),'*c');hold on%allcells(local_neighbors.columnindex',2),'MarkerSize',10,'LineWidth',1
figure(jj),text(totest(:,1)+5,totest(:,2)+5,num2str(totest(:,3)),'Color','y')
title(['At ' num2str(jj*delta_t/60) 'hrs; Track N# ' num2str(totrack) ]);

end

sametype_neighbor(jj).same = counter;
othertype_neighbor(jj).other = (nearest-counter);
fraction_same(jj).frac = counter/(nearest);
fraction_same(jj).abs = counter;
fraction_other(jj).frac = (nearest-counter)/(nearest);
fraction_other(jj).abs=(nearest-counter);
end
%close all
curr_cell_displ = cat(1,displ.dat);
curr_cell_speed = cat(1,v_2.total);
curr_cell_vx = cat(1,v.vx);
curr_cell_vy= cat(1,v.vy);
curr_cell_theta= cat(1,direction.theta);
curr_cell_rad= cat(1,direction.rad);
curr_cell_neighborhood = cat(1,fraction_same.frac);%fraction_same.frac   fraction_other.frac  
curr_cell_neighborhood2 = cat(1,fraction_other.frac);% 
same_neighbor_abs = cat(1,fraction_same.abs);
other_neighbor_abs = cat(1,fraction_other.abs);
if (toplot == 1)
% figure(4),scatter(curr_cell_speed,curr_cell_neighborhood2,[],coordintime(totrack).dat(1:end-1,time),'filled','Marker','p');box on;ylabel('Fraction of neighbors of the other cell type');xlabel('cell speed, um/hr');title('Color:time, hr');hold on;colorbar
% h4 = figure(4);h4.Colormap = jet;ylim([0 1]);%caxis([0 1]);
figure(2),polarscatter(curr_cell_rad(1:x2),curr_cell_speed(1:x2),coordintime(totrack).dat((1:x2),time),curr_cell_neighborhood(1:x2),'filled','MarkerEdgeColor','k');hold on%coordintime(totrack).dat(1:end-1,time)
figure(3),plot(coordintime(totrack).dat(1:end-1,time),curr_cell_neighborhood,'-kp','MarkerFaceColor',colormap(randi(size(colormap,1)),:),'MarkerSize',10);hold on
figure(5),scatter(curr_cell_displ(1:x2),curr_cell_neighborhood2(1:x2),[],coordintime(totrack).dat(1:x2,time),'filled','Marker','p');box on;title('Color:Time, frames');ylabel('Fraction of neighbors of the other cell type ');xlabel('Cell displacement,um');hold on;colorbar
h5 = figure(5);h5.Colormap = jet;ylim([0 1]);
cc = corrcoef(curr_cell_displ(1:x2),curr_cell_neighborhood2(1:x2));
cc = cc(1,2);
text(curr_cell_displ(end-1),0.9,['Correlation coefficient ' num2str(cc)]);
h = figure(2);
%ylim([0 1.1]); xlim([0 max(curr_cell_speed)]);%max(curr_cell_velocity)
box on
str1 = "Color: Fraction of same type neighbors within " + num2str(round(userParam.local_sz*userParam.pxtomicron))+"um neighborhood";
titlestr = "Velocity of Trackedcell wrt +X direction, r(um/hr) theta(degrees)" + "\n" + str1 + " \n"+"Label size increases with frame number";
titlestr = compose(titlestr);
title(titlestr);
h.Colormap = jet;
caxis([0 1]);
colorbar
h1 = figure(3);
ylim([0 1.1]);%xlim([0 max(curr_cell_velocity)]);
box on
xlabel('Time, hours');
ylabel('Fraction of neighbors of the same cell type')
title(['Tracked CFP cell; Neighborhood size ' num2str(round(userParam.local_sz*userParam.pxtomicron)) 'um; Track N# ' num2str(totrack)]);
h1.Colormap = jet;
X = coordintime(totrack).dat(end,time);
h1.CurrentAxes.XTick = (1:7:X);
h1.CurrentAxes.XTickLabel = (1:7:X)*delta_t/60;
end
end
samecell_neighbors=curr_cell_neighborhood;
othercell_neighbors=curr_cell_neighborhood2;


end
