% diffusion properties from msd
% 1. Get the D coeff for chuncks of time , e.g. 4 dts and average the
% fractions of unlike/like cells in those time intervals; then plot such D
% vs fractins
% 2. Fit the full msd vs time curve with the multi term equation for free
% and directed motion to ge the D and v;
clear all
close all
str1 = '200umR, pluri cells' ;
N = 3;
delta_t = 15;
str2 = ['SubDiff,200umR,pluri cells; averaging every ' num2str(N*delta_t/60) 'hrs'];

matfile1={'qMotion_f0009_w0001_time_0.25hrs-10.5hrs.mat','qMotion_f00010_w0001_time_0.25hrs-11.25hrs.mat','qMotion_f00011_w0001_time_0.25hrs-11.25hrs.mat','qMotion_f00012_w0001_time_0.25hrs-10.5hrs.mat'};%,'anom_D','confined_motion','diff_flow'
matfile2 ={'SelectedTrackIDsData_f9_chan_w0001.mat','SelectedTrackIDsData_f10_chan_w0001.mat','SelectedTrackIDsData_f11_chan_w0001.mat','SelectedTrackIDsData_f12_chan_w0001.mat'};
x3= 1;
nbr_same = 0;
choosetracktype = 4;

%   confined_motion: choosetracktype = 1;  
%   anom_D:   choosetracktype = 2;
%   diff_flow:   choosetracktype = 3;
%   msd_all:    choosetracktype = 4;%  
%   confined_motion =[]
%   anom_D =[];
%   diff_flow=[];
D_all2=struct;
nbrh_all2=struct;
time2=struct;
for i=1:size(matfile1,2)
[D_all,nbrh_all,time]=get_trackspecific_data(matfile1{i},matfile2{i},N,x3,nbr_same,choosetracktype);
D_all2(i).dat = D_all;
nbrh_all2(i).dat =  nbrh_all;
time2(i).dat =  time;
end

t1 = [D_all2.dat];
t2=  [nbrh_all2.dat];
t3=  [time2.dat];
szall=[t1.sz];
figure(2), title(str2);ylim([ 0 4000])
%%
close all
tmp3 = szall;
maxsz=max(tmp3);
timeintervals=[];
diff_in_t=[];
nbrh_in_t=[];
q = 1;
for jj=1:size(t1,2)
    if ((t1(jj).sz) == maxsz)
timeintervals(1:maxsz,1) = t3(jj).dat.*delta_t/60;        
diff_in_t(1:maxsz,q) = t1(jj).dat;
nbrh_in_t(1:maxsz,q)=t2(jj).dat;
 q = q+1;
    end   
end
% diff_in_t=[D_all.dat];
% nbrh_in_t = [nbrh_all.dat];

meanD=zeros(size(diff_in_t,1),1);
errD=zeros(size(diff_in_t,1),1);
meanNbrh=zeros(size(diff_in_t,1),1);
for jj=1:size(diff_in_t,1)
 %figure(6), scatter((timevect(isfinite(nbr_vardt))')*delta_t/60,abs(var(isfinite(nbr_vardt))),[],nbr_vardt(isfinite(nbr_vardt))','filled');box on;hold on%'p','Markersize',12); hold on;box on
 %hold on,figure(6), errorbar(timeintervals(jj),mean(abs(diff_in_t(jj,:))),std(abs(diff_in_t(jj,:))),'xk','LineWidth',1.2);hold on
meanD(jj)=mean(abs(diff_in_t(jj,:)));
errD(jj)=std(abs(diff_in_t(jj,:)));
meanNbrh(jj)=mean(nbrh_in_t(jj,:));
end
figure(8),errorbar(timeintervals,meanD,errD,'k.');hold on
scatter(timeintervals,meanD,[],meanNbrh,'filled');box on;hold on
h8 = figure(8);
h8.Colormap = jet;
h8.CurrentAxes.LineWidth = 2;
h8.CurrentAxes.FontSize = 10;
colorbar
caxis([0.2 0.7]);%(min(meanNbrh)) (max(meanNbrh))
ylim([(min(errD)-max(meanD)+5) max(errD)+max(meanD)+5])
xlim([0 timeintervals(end)+1])
xlabel('Time, hours ');
ylabel('D, um^2/hr, (each pt is <D> over tracks at that time pt');
if nbr_same == 1
title(['Color: <same type neighborhood (fraction)>; |D|, ' str1 ])
else
title(['Color: < unlike neighborhood (fraction)>;  |D|, ' str1 ])
end

%% unused for now

 paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
str1 = '200umR, CFP cells';
run(paramfile)
global userParam
colormap2 = prism;
cc = struct;
time = struct;
cc3 = struct;
%figure(1),plot(msd_all(xx).dat(1:x2,2),msd_all(xx).dat(1:x2,1),'-.','MarkerFaceCOlor',colormap2(xx,:),'Markersize',14,'LineWidth',2); box on
%close all
N = 3;
x1 =1;
%x2 =45;
nbr_same =0;
%   confined_motion =[];%  choosetracktype = 1;  
%   anom_D = [] ;%   choosetracktype = 2;
%   diff_flow =[];%  choosetracktype = 3;
%   diff_and_aD = cat(2,diff_flow,anom_D);
%    msd_all    choosetracktype = 4;
D_all = struct;
nbrh_all = struct;
choosetracktype = 1;

for xx = 1:size(diff_flow,2)%size(diff_flow,2)%size(msd_all,2)
xx =diff_flow(xx);
% for 2D trajectory, MSD = 4Dt;  get the slope from the fit D=slope/4 [um^2/hour]
ytofit_all = msd_all(xx).dat(x1:x2,1);
xtofit_all = ((x1:x2)')*delta_t/60;
nbr_other_all = allcells_local(xx).other(x1:x2);
nbr_other_all=nbr_other_all(isfinite(nbr_other_all));
nbr_same_all = allcells_local(xx).same(x1:x2);
nbr_same_all=nbr_same_all(isfinite(nbr_same_all));
displ_all= allcells_disp_inT(xx).dat(x1:x2);
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

for k=x1:N:(new_sz-N)    
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

    if k == (new_sz-N)
    text(tmp2(end)+0.1,tmp1(end),['trackID' num2str(goodtracks(xx)) '(' num2str(xx) ')' ],'color',colormap2(c,:),'FontSize',10);
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
title(str1);
figure(3), plot(nbr_vardt',var,'p','MarkerFaceColor',colormap2(c,:),'MarkerSize',10);box on; hold on
%scatter(nbr_vardt',cat(1,D_vardt.dat),[],timevect);box on; hold on;colorbar
D_all(xx).dat =var; 
D_all(xx).sz =size(var,1);
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

figure(6), scatter((timevect(isfinite(nbr_vardt))')*delta_t/60,abs(var(isfinite(nbr_vardt))),[],nbr_vardt(isfinite(nbr_vardt))','filled');box on;hold on%'p','Markersize',12); hold on;box on
h6 = figure(6);
h6.Colormap = jet;
caxis([0 1]);colorbar
ylim([ 0 max(cat(1,D_all.dat))+5])%-(max(cat(1,D_all.dat))+5)
xlabel(['Time, hours (each pt is average over)' num2str(N*delta_t/60) 'hrs']);
ylabel(['D, mu^2/hr, (each pt is average over)' num2str(N*delta_t/60) 'hrs; ' ]);
if nbr_same == 1
title(['Color: averaged same type neighborhood (fraction); |D|, '  str1])
else
title(['Color: averaged unlike neighborhood (fraction); |D|,'  str1 ])
end
end
%-------------
tmp3 = [D_all.sz];
maxsz=max(tmp3);
timeintervals=[];
diff_in_t=[];
nbrh_in_t=[];
q = 1;
for jj=1:size(D_all,2)
    if ((D_all(jj).sz) == maxsz)
timeintervals(1:maxsz,1) = time(jj).dat.*delta_t/60;        
diff_in_t(1:maxsz,q) = D_all(jj).dat;
nbrh_in_t(1:maxsz,q)=nbrh_all(jj).dat;
 q = q+1;
    end   
end
% diff_in_t=[D_all.dat];
% nbrh_in_t = [nbrh_all.dat];

meanD=zeros(size(diff_in_t,1),1);
errD=zeros(size(diff_in_t,1),1);
meanNbrh=zeros(size(diff_in_t,1),1);
for jj=1:size(diff_in_t,1)
hold on,figure(6), errorbar(timeintervals(jj),mean(abs(diff_in_t(jj,:))),std(abs(diff_in_t(jj,:))),'xk','LineWidth',1.2);hold on
meanD(jj)=mean(abs(diff_in_t(jj,:)));
errD(jj)=std(abs(diff_in_t(jj,:)));
meanNbrh(jj)=mean(nbrh_in_t(jj,:));
end
figure(8),errorbar(timeintervals,meanD,errD,'k.');hold on
scatter(timeintervals,meanD,[],meanNbrh,'filled');box on;hold on
h8 = figure(8);
h8.Colormap = jet;
h8.CurrentAxes.LineWidth = 2;
h8.CurrentAxes.FontSize = 10;
colorbar
caxis([(min(meanNbrh)) (max(meanNbrh))]);%(min(meanNbrh)) (max(meanNbrh))
ylim([0 max(errD)+max(meanD)+5])
xlabel('Time, hours ');
ylabel('D, um^2/hr, (each pt is <D> over tracks at that time pt');
if nbr_same == 1
title(['Color: averaged same type neighborhood (fraction); |D|' str1 ])
else
title(['Color: averaged unlike neighborhood (fraction);  |D|' str1 ])
end
%-------------
h3=figure(3);
xlim([0 1]);
ylim([0 max(cat(1,D_all.dat))+5]);
h3.CurrentAxes.LineWidth = 2;
h3.CurrentAxes.FontSize = 10;

a1=cat(1,nbrh_all.dat);
a2=cat(1,D_all.dat);
figure(7), plot(a1(isfinite(a1)),a2(isfinite(a1)),'p','MarkerSize',10);box on;xlim([0 1])
a = corrcoef(a1(isfinite(a1)),a2(isfinite(a1)));
title(['Corr coeff ' num2str(a(1,2)) '; Calculated from merged data, not per track']);
ylim([0 max(a2(isfinite(a1)))+5]);
xlim([0 1])
ylabel(['D, mu^2/hr, (each pt is average over)' num2str(N*delta_t/60) 'hrs' ]);
if nbr_same == 1
    xlabel(['Fraction of same type cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
else
    xlabel(['Fraction of unlike cells, averaged over each consecutive ' num2str(N*delta_t/60) 'hours']);
end
h7=figure(7);
h7.CurrentAxes.LineWidth = 2;
h7.CurrentAxes.FontSize = 10;

all_corr = cat(1,cc.dat);
figure(4),histogram(all_corr(isfinite(all_corr)),'normalization','probability');box on;ylim([0 1]);xlim([-1 1])
ylabel('Frequency');
if nbr_same == 1
xlabel('Corr.coeff. btw cell D and averaged same type neighborhood')
else
xlabel('Corr.coeff. btw cell D and averaged unlike neighborhood')
end
title(['Mean correlation (over tracks):' num2str(mean(all_corr(isfinite(all_corr)))) ' Averaging over each consecutive ' num2str(N*delta_t/60) 'hours']);
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
x2 = 25;
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
