% diffusion properties from msd
% 1. Get the D coeff for chuncks of time , e.g. 4 dts and average the
% fractions of unlike/like cells in those time intervals; then plot such D
% vs fractins
% 2. Fit the full msd vs time curve with the multi term equation for free
% and directed motion to ge the D and v;
 load('qMotion_f00011_w0000_time_0.25hrs-11.25hrs.mat');%,'anom_D','confined_motion','diff_flow'
% paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
% load('SelectedTrackIDsData_f11_chan_w0000.mat')
run(paramfile)
global userParam
colormap2 = prism;
cc = struct;
cc3 = struct;
%figure(1),plot(msd_all(xx).dat(1:x2,2),msd_all(xx).dat(1:x2,1),'-.','MarkerFaceCOlor',colormap2(xx,:),'Markersize',14,'LineWidth',2); box on
close all
N = 4;
%x2 =15;
nbr_same =1;
%   confined_motion =[];%    
%   anom_D = [] ;%   
%   diff_flow =[];% 
D_all = struct;
nbrh_all = struct;

for xx = 1:size(msd_all,2)%size(diff_flow,2)%size(msd_all,2)
%xx = diff_flow(xx);
% for 2D trajectory, MSD = 4Dt;  get the slope from the fit D=slope/4 [um^2/hour]
ytofit_all = msd_all(xx).dat(1:x2,1);
xtofit_all = ((1:x2)')*delta_t/60;
nbr_other_all = allcells_local(xx).other(1:x2);
nbr_other_all=nbr_other_all(isfinite(nbr_other_all));
nbr_same_all = allcells_local(xx).same(1:x2);
nbr_same_all=nbr_same_all(isfinite(nbr_same_all));
displ_all= allcells_disp_inT(xx).dat(1:x2);
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

for k=1:N:(new_sz-N)    
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
    if k == (new_sz-N)
    text(tmp2(end)+0.1,tmp1(end),['trackID' num2str(goodtracks(xx)) '(' num2str(xx) ')' ],'color',colormap2(c,:),'FontSize',10);
    end
    slope_estimate(q).slope =  cfit.p1;
    D_vardt(q).dat = slope_estimate(q).slope/4;
    disp([k k+N-1]);
    timevect(q) = k;
    q = q+1;
end
h2 = figure(2);
h2.CurrentAxes.LineWidth = 2;
h2.CurrentAxes.FontSize = 10;
ylabel('Mean Squared Displacement, um^2');
xlabel('Time, hrs');
var = cat(1,D_vardt.dat);

figure(3), plot(nbr_vardt',var,'p','MarkerFaceColor',colormap2(c,:),'MarkerSize',10);box on; hold on
%scatter(nbr_vardt',cat(1,D_vardt.dat),[],timevect);box on; hold on;colorbar
D_all(xx).dat =var; 
nbrh_all(xx).dat = nbr_vardt';

ylabel(['D coeff, calculated over each consecutive' num2str(N*delta_t/60) 'hours, um^2/hr' ]);
if nbr_same == 1
    xlabel(['Fraction of same type cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
else
    xlabel(['Fraction of unlike cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
end
cc1 = corrcoef(nbr_vardt(isfinite(nbr_vardt))',var(isfinite(nbr_vardt)));
cc(xx).dat = cc1(1,2);
title([' Neighborhood size:' num2str((userParam.local_sz*userParam.pxtomicron)) 'um']);
figure(5),plot(nbr_vardt',cell_aver_displ','p','MarkerFaceColor',colormap2(c,:),'MarkerSize',10);box on; hold on
h5=figure(5);
h5.CurrentAxes.LineWidth = 2;
h5.CurrentAxes.FontSize = 10;
ylabel(['Cell displacement, averaged over each consecutive' num2str(N*delta_t/60) 'hours, um^2/hr' ]);
if nbr_same == 1
    xlabel(['Fraction of same type cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
else
    xlabel(['Fraction of unlike cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
end
cc2 = corrcoef(nbr_vardt(isfinite(nbr_vardt))',cell_aver_displ(isfinite(nbr_vardt))');
cc3(xx).dat = cc2(1,2);

figure(6), scatter((timevect(isfinite(nbr_vardt))')*delta_t/60,(var(isfinite(nbr_vardt))),[],nbr_vardt(isfinite(nbr_vardt)),'filled');box on;hold on%'p','Markersize',12); hold on;box on
caxis([0 1]);colorbar
xlabel(['Time, hours (each pt is average over)' num2str(N*delta_t/60) 'hrs']);
ylabel(['D, mu^2/hr, (each pt is average over)' num2str(N*delta_t/60) 'hrs' ]);
if nbr_same == 1
title('Color: averaged same type neighborhood')
else
title('Color: averaged unlike neighborhood')
end
end


h3=figure(3);
xlim([0 1]);
ylim([0 max(cat(1,D_all.dat))+5]);
h3.CurrentAxes.LineWidth = 2;
h3.CurrentAxes.FontSize = 10;

a1=cat(1,nbrh_all.dat);
a2=cat(1,D_all.dat);
figure(7), plot(a1(isfinite(a1)),a2(isfinite(a1)),'p','MarkerSize',10);box on;xlim([0 1])
a = corrcoef(a1(isfinite(a1)),a2(isfinite(a1)));
title(['Corr coeff ' num2str(a(1,2)) '; All data merged ']);
%'Correlation coefficient:' num2str(cc(xx).dat)
all_corr = cat(1,cc.dat);
figure(4),histogram(all_corr(isfinite(all_corr)),'normalization','probability');box on;ylim([0 1]);xlim([-1 1])
ylabel('Frequency');
if nbr_same == 1
xlabel('Corr.coeff. btw cell D and averaged same type neighborhood')
else
xlabel('Corr.coeff. btw cell D and averaged unlike neighborhood')
end
title(['Mean correlation:' num2str(mean(all_corr(isfinite(all_corr)))) ' Averaging over each consecutive ' num2str(N*delta_t/60) 'hours']);
figure(5)
all_corr2 = cat(1,cc3.dat);
title([' Neighborhood size:' num2str((userParam.local_sz*userParam.pxtomicron)) 'um Corr.: ' num2str(mean(all_corr2(isfinite(all_corr2))))]);

%% fit the full msd vs time track to model
% here a is the diffusion coefficient D
% 4*a*power(x,alpha)  alpha<1, anomalous diffusion
% 4*a*x + (v*x)*(v*x) directed motion with iffusion
% 4*a*t normal diffusion
close all
% fit diffusion+flow
D_fromfit=[];
v_fromfit=[];
q = 1;
x2 = 45;
for  xx = 1:size(diff_flow,2)%size(msd_all,2) confined_motion diff_flow
    xx = diff_flow(xx);   
    ytofit_all = msd_all(xx).dat(1:x2,1);
    x = ((1:x2)')*delta_t/60;
    myfittype =  fittype('4*a*x + (v*x)*(v*x)','dependent',{'ytofit_all'},...
        'independent',{'x'},'coefficients',{'a','v'}) ;
    myfit = fit(x,ytofit_all,myfittype,'StartPoint',[20 10]);
    figure(4),plot(x,ytofit_all,'mp');hold on
    plot(myfit,'k');hold on
    ylabel('MSD, um^2')
    xlabel('time, hrs')
    legend('off')
    D_fromfit(q,1)=myfit.a;
    v_fromfit(q,1)=myfit.v;
    q=q+1;
end
figure, histogram(D_fromfit,'normalization','probability');title('D, um^2/hr');ylim([0 1])
figure, histogram(v_fromfit,'normalization','probability');title('v, um/hr');ylim([0 1])
% figure,scatter(1:size(v_fromfit,1),v_fromfit)
%% fit anomalous diffusion
D_fromfit=[];
alpha_fromfit=[];
close all
q = 1;
x2 = 45;
for xx = 1:size(anom_D,2)%size(msd_all,2)anom_D
    xx =anom_D(xx);    
    ytofit_all = msd_all(xx).dat(1:x2,1);
    x = ((1:x2)')*delta_t/60;
    myfittype =  fittype('4*a*power(x,alpha)','dependent',{'ytofit_all'},...
        'independent',{'x'},'coefficients',{'a','alpha'}) ;
    myfit = fit(x,ytofit_all,myfittype,'StartPoint',[20 0.5]);
    figure(4),plot(x,ytofit_all,'mp');hold on
    plot(myfit,'k');hold on
    ylabel('MSD, um^2')
    xlabel('time, hrs')
    legend('off')
    D_fromfit(q,1)=myfit.a;
    alpha_fromfit(q,1)=myfit.alpha;
    q=q+1;
end
figure, histogram(D_fromfit,'normalization','probability');title('D, um^2/hr');ylim([0 1])
figure, histogram(alpha_fromfit,'normalization','probability');title('alpha, um/hr');ylim([0 1])
%% anomaluos diffusion with flow
% TODO: restrict D to be positive
% msd = msd0*(1-A*exp(-4*B*a*t/msd0))
D_fromfit=[];
alpha_fromfit=[];
v_fromfit=[];
clear myfit
close all
q = 1;
x2 = 45;
% StartPOint option in fit: initial valuoes for the coefficients
% 
for xx  =1%:size(diff_flow,2)%size(msd_all,2)diff_flow  anom_D
    xx = diff_flow(xx);     
    ytofit_all = msd_all(xx).dat(1:x2,1);
    x = ((1:x2)')*delta_t/60;
    myfittype =  fittype('4*a*power(x,alpha)+power(v*x,2)','dependent',{'ytofit_all'},...
        'independent',{'x'},'coefficients',{'a','alpha','v'}) ;
    myfit = fit(x,ytofit_all,myfittype,'StartPoint',[20 0.5 10]);
    figure(4),plot(x,ytofit_all,'mp');hold on
    plot(myfit,'k');hold on
    ylabel('MSD, um^2')
    xlabel('time, hrs')
    legend('off')
    D_fromfit(q,1)=myfit.a;
    alpha_fromfit(q,1)=myfit.alpha;
    v_fromfit(q,1)=myfit.v;
    q=q+1;
end
figure, histogram(D_fromfit,'normalization','probability');title('D, um^2/hr');ylim([0 1])
figure, histogram(alpha_fromfit,'normalization','probability');title('alpha, um/hr');ylim([0 1])
figure, histogram(v_fromfit,'normalization','probability');title('v, um/hr');ylim([0 1])
figure,scatter(1:size(alpha_fromfit,1),alpha_fromfit,'p','filled');title('alpha, um/hr');box on

%% fit free diffusion
D_fromfit=[];
alpha_fromfit=[];
close all
q = 1;
x2 = 40;
for xx = 2%:size(anom_D,2)%size(msd_all,2)
    xx = confined_motion(xx);    
    ytofit_all = msd_all(xx).dat(1:x2,1);
    x = ((1:x2)')*delta_t/60;
    myfittype =  fittype('4*a*x','dependent',{'ytofit_all'},...
        'independent',{'x'},'coefficients',{'a'}) ;
    myfit = fit(x,ytofit_all,myfittype,'StartPoint',30);
    figure(4),plot(x,ytofit_all,'mp');hold on
    plot(myfit,'k');hold on
    ylabel('MSD, um^2')
    xlabel('time, hrs')
    legend('off')
    D_fromfit(q,1)=myfit.a;
    %alpha_fromfit(q,1)=myfit.alpha;
    q=q+1;
end
