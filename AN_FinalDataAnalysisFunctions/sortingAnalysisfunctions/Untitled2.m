%% load the raw images data for selected channel
pos = 17;
chan = 1;
imagedir  ='C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2017-07-14-Smad4sorting_maxProjections';
ff1 = readAndorDirectory(imagedir); 
nucmoviefile = getAndorFileName(ff1,ff1.p(pos),[],[],ff1.w(chan));%getAndorFileName(files,pos,time,z,w)
% nreader = bfGetReader(nucmoviefile);
% nt = nreader.getSizeT;
r3 = bfopen(nucmoviefile);

%% find numeric ID of the long tracks
q = 1;
longtracks = [];
for jj=1:size(coordintime,2)
tmp = size(coordintime(jj).dat,1);
if tmp >=25
longtracks(q,1) = jj;
end
q = q+1;
end
goodtracks_tmp = nonzeros(longtracks);

%% make movie of the track for the track number 'totrack' to check that the track corresponds to actual cell and good
celltype = 'CFP';%CFP pluri
matlabtracking = 1;%
%pos14 [9,11,19,21,30,36,40,41,57,59*,60,64,65,67,86,88,93,94,100,102,104,133,144,145,147,148,149]
%pos16_w1: [3,13,20,27,36,40*,41*,45,51*,52*,55,60*,63*]
totrack = 67;
%good track IDs pos0chanPluri:13,
%,36,40,45,49,52,64,79,84,85,93,96,104,108,110,111*,129,131,143,153,160,171,174,182,183,196,199,204,214,224*,230,
% 231,232,
%good track IDs pos0chanCFP:
%[17,19,20,24,33,34,38,44,50,53,56,57,58,61,68,69,72,75,83,88,89,104,105,110,125,128]
nt = 80; % segmented only 80 tpts
trackcelltype_movie = ViewTrack_getMovie(coordintime,celltype,matlabtracking,totrack,r3);
%% run the movie
close all
h = figure(1);%h.Colormap = jet; %
movie(h,trackcelltype_movie,1,3);%,[0 0 560 420]
% save the movie
% pos1 = ff1.p(pos);
% pos_video = VideoWriter(['C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\Ilastik_AutoTracked_movies\TrackSortingCells_Movie_pos_' num2str(pos1) ' Track# ' num2str(totrack) ' CellType_' celltype '.avi']);%C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\For_MatlabTracking\Sorting_Movie_with_CellTrack
% pos_video.FrameRate = 2;
% open(pos_video);
% writeVideo(pos_video,trackcelltype_movie);
% close(pos_video);
%%