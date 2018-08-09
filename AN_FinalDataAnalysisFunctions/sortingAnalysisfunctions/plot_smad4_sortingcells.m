% plot the data from GUI (Smad4 nuclear-to-cyto ratio for selected cells in time)
close all
ff = dir('.');
delta_t = 15;
title_str = 'Inside cells only';
%colSZ = 200 um
for ii =1:size(ff,1)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'SelectedCellsData'))
matfile = ff(ii).name;
load(matfile);
disp(matfile)
time_vect = ((fluor_selectedtracks.tp).*delta_t)/60;
s4_vect = fluor_selectedtracks.fluor;
figure(1),plot(time_vect, s4_vect,'bp','Markersize',9); hold on
box on
    end
    
end
figure(1),plot(unique(time_vect),ones(size(unique(time_vect))),'--r','LineWidth',3);
h = figure(1);
h.CurrentAxes.LineWidth = 2;
h.CurrentAxes.FontSize = 12;
ylabel('Smad4 nuclear to cyto ratio, a.u.')
xlabel('Time,hours')
title(['Cells picked randomly, not tracked data. ' title_str ]);
ylim([0.3 1.7])