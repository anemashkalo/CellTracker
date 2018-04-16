
function [D_all,nbrh_all,time]=get_trackspecific_data(matfile1,matfile2,N,x3,nbr_same,choosetracktype)

%function to extract the data from the 'qMotion..' file for each type of msd tracks  

%   confined_motion =[];%  choosetracktype = 1;  
%   anom_D = [] ;%   choosetracktype = 2;
%   diff_flow =[];%  choosetracktype = 3;
%   diff_and_aD = cat(2,diff_flow,anom_D);
%    msd_all    choosetracktype = 4;
%   confined_motion =[];%    
%   anom_D = [] ;%   
%   diff_flow =[];
%    msd_all
clear msd_data
load(matfile1);
load(matfile2);
disp(['Loaded : matfile:' matfile1 '; and : ' matfile2 ])
D_all = struct;
nbrh_all = struct;
colormap2 = jet;
if choosetracktype == 1
    msd_data = confined_motion;
    
elseif choosetracktype == 2
    msd_data = anom_D;
    
elseif choosetracktype == 3
    msd_data = diff_flow  ;
    
elseif choosetracktype == 4
    msd_data = msd_all;
    
end
            
for xx = 1:size(msd_data,2)%size(diff_flow,2)%size(msd_all,2)
   if ~isempty(msd_data)
    if choosetracktype < 4 % 
ii =msd_data(xx);
    end
    if choosetracktype == 4 % 
ii =xx;%goodtracks(xx)
    end
    % for 2D trajectory, MSD = 4Dt;  get the slope from the fit D=slope/4 [um^2/hour]
ytofit_all = msd_all(ii).dat(x3:x2,1);
xtofit_all = ((x3:x2)')*delta_t/60;
nbr_other_all = allcells_local(ii).other(x3:x2);
nbr_other_all=nbr_other_all(isfinite(nbr_other_all));
nbr_same_all = allcells_local(ii).same(x3:x2);
nbr_same_all=nbr_same_all(isfinite(nbr_same_all));
displ_all= allcells_disp_inT(ii).dat(x3:x2);
tmp1=[];
tmp2 = [];
cell_aver_displ=[];
fit_out=struct;
slope_estimate=struct;
D_vardt=struct;
timevect=[];
q = 1;
nbr_vardt=[];
c = randi(size(colormap2,1));
if nbr_same == 1
    new_sz = size(nbr_same_all,1);
else
   new_sz = size(nbr_other_all,1) ;
end

for k=x3:N:(new_sz-N)    
    if nbr_same == 1
    ytofit_all = ytofit_all(isfinite(nbr_same_all));
    xtofit_all = xtofit_all(isfinite(nbr_same_all));
    tmp1 = ytofit_all(k:k+N-1);
    tmp2 = xtofit_all(k:k+N-1);    
    nbr_vardt(q)=mean(nbr_same_all(k:k+N-1));%nbr_other_all
    else
    ytofit_all = ytofit_all(isfinite(nbr_other_all));
    xtofit_all = xtofit_all(isfinite(nbr_other_all));
    tmp1 = ytofit_all(k:k+N-1);
    tmp2 = xtofit_all(k:k+N-1);       
    nbr_vardt(q)=mean(nbr_other_all(k:k+N-1));%nbr_other_all
    end
    cell_aver_displ(q) = mean(displ_all((k:k+N-1)));
    cfit = fit(tmp2,tmp1,'poly1');
    fit_out(q).dat = cfit;
    figure(2),plot(tmp2,tmp1,'p','MarkerSize',8,'Color',colormap2(c,:));hold on
    plot(tmp2,cfit(tmp2),'LineWidth',1.5,'Color','k');hold on
      %  plot(cfit);hold on

    if (k == (new_sz-N-2)) || (k == (new_sz-N-1))|| (k == (new_sz-N))
        if choosetracktype < 4
    text(tmp2(end)+0.1,tmp1(end),['trackID' num2str(msd_data(xx)) '(' num2str(xx) ')' ],'color',colormap2(c,:),'FontSize',10);
        else
    text(tmp2(end)+0.1,tmp1(end),['trackID' num2str(goodtracks(xx)) '(' num2str(xx) ')' ],'color',colormap2(c,:),'FontSize',10);
        end
        end
    slope_estimate(q).slope =  cfit.p1;
    D_vardt(q).dat = slope_estimate(q).slope/4;
    %disp([k k+N-1]);
    timevect(q) = k;
    time(xx).dat = timevect;
    q = q+1;
end
h2 = figure(2);
h2.CurrentAxes.LineWidth = 2;
h2.CurrentAxes.FontSize = 10;
ylabel('Mean Squared Displacement, um^2');
xlabel('Time, hrs');
var = cat(1,D_vardt.dat);
%title(str1);
% figure(3), plot(nbr_vardt',var,'p','MarkerFaceColor',colormap2(c,:),'MarkerSize',10);box on; hold on
%scatter(nbr_vardt',cat(1,D_vardt.dat),[],timevect);box on; hold on;colorbar
D_all(xx).dat =var; 
D_all(xx).sz =size(var,1);
nbrh_all(xx).dat = nbr_vardt';


   end
end

