% make projections from tif images (LS confocal slices)
% Ecad Ncad data
 chanprojection_img = [];
 allchannels_img = [];%zeros(:,:,4);
 close all
 chan = {'DAPI','RFP','CY5','CFP'};
 condition='24hr_mix_2_';
 direc2save = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-02-CAMdyn_Kong_EcadCY5_NcadRFP_CPFundiff_DAPI_RFP_CY5_CFP\maxProjections';
nz = 6;
for jj=1:size(chan,2)
    tmp = [];
for ii = 1:nz
        if ii <10
        img_now=imread(['s_C00' num2str(jj) 'Z00' num2str(ii) '.tif' ]);
        else
        img_now=imread(['s_C00' num2str(jj) 'Z0' num2str(ii) '.tif' ]);
        end
        %figure(ii),imshow(img_now,[]);
        if ii == 1
            max_img =  img_now;
        else
            max_img = max(max_img,img_now);
        end
        chanprojection_img = max_img;
        
end
%close all
 allchannels_img(:,:,jj)=chanprojection_img;
 tmp = allchannels_img(:,:,jj);
 tmp = tmp-min(tmp(:));
 tmp = tmp/max(tmp(:));
 figure(jj), imshow(tmp,[]);
 
 imwrite(tmp,[direc2save '\' condition chan{jj} '.tif']);%,'writemode','append

end
  

  
  
  