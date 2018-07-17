function [nearcellD,actual_times,dx_nearcell,dy_nearcell,dx_trackcell,dy_trackcell]=get_cell_displ(cell2,paramfile,delta_t,tracked_coord)
run(paramfile);global userParam
Displ_nearcell =struct;
dx_nearcell=struct;
dy_nearcell=struct;
dx_trackcell=struct;
dy_trackcell=struct;
nearcellD  = [];
actual_times = size(nonzeros(cell2(:,1)),1);
for jj=1:actual_times-1
        tp1 =jj;
        tp2 = jj+1;
        dt = (tp2-tp1)*delta_t/60;
        dx = (cell2(tp2,1)-cell2(tp1,1))*userParam.pxtomicron;
        dy = (cell2(tp2,2)-cell2(tp1,2))*userParam.pxtomicron;
        Displ_nearcell(jj).dat = power((dx*dx+dy*dy),0.5);
        dx_nearcell(jj).dat = dx;
        dy_nearcell(jj).dat = dy;
        dx_tr = (tracked_coord(tp2,1)-tracked_coord(tp1,1))*userParam.pxtomicron;
        dy_tr = (tracked_coord(tp2,2)-tracked_coord(tp1,2))*userParam.pxtomicron;
        dx_trackcell(jj).dat=dx_tr;
        dy_trackcell(jj).dat=dy_tr;
        
end
nearcellD = cat(1,Displ_nearcell.dat);
end