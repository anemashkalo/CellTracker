function [d_btw]=get_dist_btw_cellcenters(cell1,cell2,paramfile,delta_t)
run(paramfile);global userParam
adj_size = size(nonzeros(cell1(:,1)),1);
d_btw = zeros(adj_size,2);
for jj=1:adj_size       
    dx = (cell1(jj,1)-cell2(jj,1))*userParam.pxtomicron;
    dy = (cell1(jj,2)-cell2(jj,2))*userParam.pxtomicron;
    d_btw(jj,1)=power((dx*dx+dy*dy),0.5);
    d_btw(jj,2)=(jj*delta_t)/60;
end
end