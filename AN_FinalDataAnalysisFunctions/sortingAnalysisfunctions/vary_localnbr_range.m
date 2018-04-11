function [samecell_neighbors,othercell_neighbors,curr_cell_speed,curr_cell_displ]= vary_localnbr_range(trackID,paramfile,coordintime,delta_t,pluristats,diffstats,img_untracked,img_tracked,x2,var_sz)
run(paramfile)
global userParam
local_neighbors = [];
for totrack =trackID
tpts = size(coordintime(totrack).dat,1);%size(cfpstats,2)  
v_2 = struct; % total velocity at each time point
v = struct; %
othertype_neighbor = struct;
sametype_neighbor = struct;
fraction_same = struct;
fraction_other = struct;
displ = struct;
for jj = 1: (x2-1)% time points    (size(cfpstats,2)-1)
tmp1 = round(cat(1,pluristats(jj).stats.Centroid));% centrois of all cells1 at tp jj
area1 = cat(1,pluristats(jj).stats.Area);
allcells_type1 = cat(2,tmp1,area1);
tmp2 = round(cat(1,diffstats(jj).stats.Centroid));% centrois of all  cells2 at tp jj;
area2 = cat(1,diffstats(jj).stats.Area);
allcells_type2 = cat(2,tmp2,area2);
curr_track = coordintime(totrack).dat(jj,:);% coord of cell of other type, for which the neighborhood is quantified
rawimg1 = (img_untracked{1}{jj});% 
rawimg1 = imadjust(rawimg1,stretchlim(rawimg1));
rawimg = (img_tracked{1}{jj});% 
rawimg = imadjust(rawimg,stretchlim(rawimg));
total_img = max(rawimg1,rawimg);
%at each time point find how many cells of each type are surrounding
% the given cell ( within given distance in pixels)
allcells = cat(1,allcells_type1,allcells_type2); % all cells, including the coord of a current cell 
local_neighbors = ipdm(coordintime(totrack).dat(jj,1:2),allcells(:,1:2),'Result','Structure','Subset','Maximum','Limit',var_sz);% get all the cells, closer than local_sz
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
% if (jj == 1) %last tp: (tpts-1) or x2
% %figure(jj),imshow(total_img,[500 1500]);hold on % 
% figure(jj),imshowpair(rawimg1,rawimg);hold on % show both cell types on one image
% figure(jj),plot(coordintime(totrack).dat(jj,1),coordintime(totrack).dat(jj,2),'kp','MarkerFaceColor','y','MarkerSize',12,'LineWidth',1);hold on%colormap(randcolor,:)
% figure(jj),plot(totest(:,1),totest(:,2),'*c');hold on%allcells(local_neighbors.columnindex',2),'MarkerSize',10,'LineWidth',1
% figure(jj),text(totest(:,1)+5,totest(:,2)+5,num2str(totest(:,3)),'Color','y')
% title(['At ' num2str(jj*delta_t/60) 'hrs; Track N# ' num2str(totrack) ]);
% end
sametype_neighbor(jj).same = counter;
othertype_neighbor(jj).other = (nearest-counter);
fraction_same(jj).frac = counter/(nearest);
fraction_same(jj).abs = counter;
fraction_other(jj).frac = (nearest-counter)/(nearest);
fraction_other(jj).abs=(nearest-counter);
end
%close all
% curr_cell_vx = cat(1,v.vx);
% curr_cell_vy= cat(1,v.vy);
end
curr_cell_displ = cat(1,displ.dat);
curr_cell_speed = cat(1,v_2.total);
samecell_neighbors=cat(1,fraction_same.frac);%fraction_same.frac   fraction_other.frac ;
othercell_neighbors=cat(1,fraction_other.frac);% 

end