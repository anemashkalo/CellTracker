function [sametype,othertype,curr_track,Displ_trackedcell]=get_local_cell_stats(coordintime,trackID,untrackedstats,trackedstats,x2,paramfile,delta_t,toplot,toshow)
sametype = struct;
othertype= struct;     
Displ_trackedcell = struct;
run(paramfile)
global userParam
for totrack =trackID    
    for jj = 1: (x2-1)% time points    (size(cfpstats,2)-1)
        % disp(jj)
        coordinate_same = struct;
        coordinate_other = struct;
        coordinate_same = [];
        coordinate_other=[];
        tmp1 = round(cat(1,untrackedstats(jj).stats.Centroid));% centrois of all pluri cells at tp jj
        area1 = cat(1,untrackedstats(jj).stats.Area);
        allcells_type1 = cat(2,tmp1,area1);
        tmp2 = round(cat(1,trackedstats(jj).stats.Centroid));% centrois of all pluri cells at tp jj;
        area2 = cat(1,trackedstats(jj).stats.Area);
        allcells_type2 = cat(2,tmp2,area2);
        
        allcells = cat(1,allcells_type1,allcells_type2); % all cells, including the coord of a current cell
        local_neighbors = ipdm(coordintime(totrack).dat(jj,1:2),allcells(:,1:2),'Result','Structure','Subset','Maximum','Limit',userParam.local_sz);% get all the cells, closer than local_sz
        tp1 =jj;
        tp2 = jj+1;
        % disp([tp1 tp2]);
        dt = (tp2-tp1)*delta_t/60;
        dx = (coordintime(totrack).dat(tp2,1)-coordintime(totrack).dat(tp1,1))*userParam.pxtomicron;
        dy = (coordintime(totrack).dat(tp2,2)-coordintime(totrack).dat(tp1,2))*userParam.pxtomicron;
        curr_track(jj,1:3) = coordintime(totrack).dat(jj,1:3);
        Displ_trackedcell(jj).dat = power(dx*dx+dy*dy,0.5);
        
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
        q1 = 1;
        q2 = 1;
        
        for h=1:size(totest,1)
            % if one of the found local neighborhood cells are in the set of
            % celltype2, the nearest neighbor to it will be at zero distance (the cell is closest to itself);
            % counting how many of those, gives the number of neighbors of cell
            % type2, the rest (nearest-counter) is the other cell type
            tmp = ipdm(totest(h,1:2),allcells_type2(:,1:2),'Result','Structure','Subset','NearestNeighbor');%
            if tmp.distance == 0
                %disp('same')
                coordinate_same(q1).xya = totest(h,1:3);% save the coordinate and area of the nearest like cell
                q1 = q1+1;
            else
                %disp('other')
                coordinate_other(q2).xya = totest(h,1:3);% save the coordinate and area of the nearest unlike cell
                q2 = q2+1;
            end
        end
        
        if size(coordinate_same,2) == size(totest,1)% if all cells are same type
            coordinate_other = [];
        end
        
        if size(coordinate_other,2) == size(totest,1)% if all cells are other type
            coordinate_same = [];
        end
        if toplot == 1 && jj == toshow
            if ~isempty(coordinate_same)
                tmp_same = cat(1,coordinate_same.xya);
            end
            if ~isempty(coordinate_other)
                tmp_other = cat(1,coordinate_other.xya);
            end
            figure(1), imshowpair(untrackedstats(jj).img,trackedstats(jj).img);hold on
            plot(coordintime(totrack).dat(jj,1),coordintime(totrack).dat(jj,2),'pc','Markersize',8,'MarkerFaceColor','c','MarkerEdgeColor','r');hold on
            if ~isempty(coordinate_same)
                plot(tmp_same(:,1),tmp_same(:,2),'co','Markersize',8,'MarkerFaceColor','c');hold on
            end
            if ~isempty(coordinate_other)
                plot(tmp_other(:,1),tmp_other(:,2),'bo','Markersize',8,'MarkerFaceColor','b');hold on%,'LineWidth',0.8
            end
            title(['Time point  ' num2str(delta_t*jj/60) ' hrs; frame: ' num2str(jj) ]);
        end
        if ~isempty(coordinate_same)
            sametype(jj).dat = cat(1,coordinate_same.xya);
        else
            sametype(jj).dat=[];
        end
        if ~isempty(coordinate_other)
            othertype(jj).dat  = cat(1,coordinate_other.xya);
        else
            othertype(jj).dat=[];
        end
    end
end

end