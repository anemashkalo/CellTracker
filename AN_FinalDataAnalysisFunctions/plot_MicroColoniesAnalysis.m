%% plot mean valuoes of protein expression for each experimental condition
 close all
 nms = {'esi_GNANOG_RSOX2','nodal_ko_GNANOG_RSOX2'}; % DAPI GFP(nanog) RFP(sox2) 
 nms2 = {'esi17 cells','Nodal knockout cells'};
 dapimax = 60000;% if want to look at cells only below the threshold of dapi intensity
 chanmax = 60000; %if want to look at cells only below the threshold of dapi intensity
 dir = '.';
 %------------ rescale the dapi mean values within the dataset
 scaledapi =1; % set to 1 if want to rescale, set to 0 if not.  
 index1 = [5]; % dapi channel, do not change
 if scaledapi == 1
for k=1:size(nms2,2)
[dapi(k),ncells] = getmeandapi(nms(k),dir,index1, dapimax);
disp(['cells found' num2str(ncells) ]);
disp(['Mean dapi' num2str(dapi(k)) ]);
end
dapiscalefactor = dapi/dapi(1);
end
if scaledapi == 0 % 
dapiscalefactor = ones(1,size(nms,2));
end
%----------------
 indexvar = [6 8]; % channels in peaks structure that you want to look at mean values of
 param1 = {'Nanog','Sox2'};% what the channels represent
 colormap = colorcube; 
 C = {'b','m','g','r'};    
 chans = size(indexvar,2);
 vect = [1:size(nms2,2)];
 newdata = cell(1,size(param1,2)); % will be returned below, contains mean and standard deviation
 %GeneralizedMeanAN(nms,nms2,dir,midcoord,fincoord,index1,param1,plottype,flag,dapimax, chanmax,dapiscalefactor)
 for k=1:chans
 [newdata{k}] = GeneralizedMeanAN(nms,nms2,dir,[],[],[indexvar(k) 5],param1{k},0,0,dapimax,chanmax,dapiscalefactor);
   figure(1),errorbar(vect',newdata{k}(:,1),newdata{k}(:,2),'-.','color',C{k},'markersize',16,'linewidth',1.3);hold on
   max_y(k) = max(newdata{k}(:,1)); 
   min_y(k) = min(newdata{k}(:,1)); 
  end
 title('Community effect in microcolonies') ; hh = figure(1);box on
 hh.CurrentAxes.LineWidth = 3; hh.CurrentAxes.FontSize = 12; ylabel('mean/dapi');
 hh.CurrentAxes.XTick = vect;
 hh.CurrentAxes.XTickLabel =  nms2 ;
 hh.CurrentAxes.XTickLabelRotation = -15;
 hh.CurrentAxes.XLim = [vect(1)-1 vect(end)+1];
 hh.CurrentAxes.YLim = [min(min_y)-0.5 max(max_y)+0.5];
% hh.CurrentAxes.Color = colormap; 
 legend(param1,'Orientation','horizontal');
 %% plot colony size dependencies    (this section generates 4 figures)
 % run for each channel independently (specify index1, param1, thresh, since those are not the same for each channel)
 
 index1 = [8 5]; % which data from peaks structure in the matfile to use, index1(1) = stain; index1(2) = dapi, always 5;
 param1 = 'Sox2';% what protein expression or other this channel represents (corresponds to index1(1))
 thresh =[2.3 2.3];   % what is considered positive expression in the channel currently analyzed (normalized to dapi) 
 % this value needs to be estimated from images of the positive control,
 % for example (take a improfile from a cell in dapi and channel, divide, this will be an estimate for the thresh)
 % you can also take a mean value for the channel, generated in the script
 % section above 
 
 flag1 = 1;      % if flag1 and flag are = 0, will not generate plots, just return calculated variables
 flag = 1; 
  [totalcells,ratios,ratios2,totcol] = GeneralizedColonyAnalysisAN(thresh,dir,nms,nms2,[],[],index1,param1,0,flag,dapimax,chanmax,dapiscalefactor);
  [rawdata] =  Intensity_vs_ColSize(nms,nms2,dir,index1,param1,dapimax,chanmax,dapiscalefactor,flag1);
 
%% plot the scatter plots colorcoded by colonysize
param1 = 'Sox2'; % DAPI GFP(nanog) RFP(sox2)
param2 = 'Nanog';
index1 = [5];
toplot = cell(1,size(nms,2));
flag = 0;% generate third column with the colony size
flag2 = 1;% do not normalize to DAPI if flag == 0;
index2 = [8 6]; % DAPI GFP(nanog) RFP(sox2)
ucol = 5;% max colony size that want to look at
for k=1:size(nms,2)
        filename{k} = [dir filesep  nms{k} '.mat'];
        load(filename{k},'peaks','dims','plate1');
        col = plate1.colonies;
[alldata] = mkVectorsForScatterAN(peaks,col,index2,flag,flag2,dapimax,dapiscalefactor(k),ucol);
%mkVectorsForScatterAN(peaks,col,index2,flag,flag2,dapimax,dapiscalefactor)
 toplot{k} = alldata;
end
for j=1:size(nms,2)
    figure(6+j),scatter(toplot{j}(:,2),toplot{j}(:,1),[],toplot{j}(:,3),'LineWidth',2);legend(nms2{j}); hold on %[],toplot{j}(:,3)    
    box on
    ylabel(param1)
    xlabel(param2)
    ylim([0 10]);
    xlim([0 10]);
    h=figure(6+j);
    h.Colormap = jet;
    colorbar
    caxis([1 ucol])
    title('Colorbar: colony size')
end


 