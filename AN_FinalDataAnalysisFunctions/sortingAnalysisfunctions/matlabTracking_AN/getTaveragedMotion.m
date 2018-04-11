matfile_1 = 'qMotion_f0000_w0000_time_0.25hrs-11.25hrs.mat';
str1 = 'CFP cells';
str1 = 'PLURI cells';

load(matfile_1)
N =8;
close all
cc2 = [];
for trackN = 1:size(allcells_disp_inT,2)
d = allcells_disp_inT(trackN).dat;%allcells_disp_inT
nbr_other = allcells_local(trackN).other;
tmp1=[];
tmp2 = [];
q = 1;
for k=1:N:(x2-N)
    tmp1(q) = mean(d(k:k+N-1));
    tmp2(q) = mean(nbr_other(k:k+N-1));
    disp([k k+N-1]);   
    q = q+1;
end
    cc1 = corrcoef(tmp1,tmp2);
    cc2(trackN) = cc1(2);    
    figure(1), scatter(tmp1',tmp2',[],(1:N:(x2-N+1))');box on; hold on
    h = figure(1);
    h.CurrentAxes.YLim = [0 1];
    h.Colormap = jet;
    caxis([0 x2]);
    colorbar 
    xlabel(['Cell displacement, averaged every ' num2str(N*delta_t/60) 'hrs ,um'])
    ylabel('Fraction of other cell type within 35um,same averaging')
    title([ str1 'cells; Color:Time in frames'])   

end
figure(3),histogram(cc2(isfinite(cc2)),'binwidth',0.4,'normalization','probability');box on
ylabel('Frequency');
xlabel('Corr.coeff. btw averaged cell displacement and unlike neighbothood')
ylim([0 1])
xlim([-1 1])
title(['Mean correlation coefficient:' num2str(mean(cc2(isfinite(cc2))))])
    
    