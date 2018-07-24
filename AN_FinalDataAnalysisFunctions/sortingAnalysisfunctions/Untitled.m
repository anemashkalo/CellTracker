% roughly asssociate the membrane of the cell to a cell
maxproj_dir = 'C:\Users\Nastya\Desktop\RiceResearch\2017-10-04-REMOTE_WORK\2018-07-02-CAMdyn_Kong_EcadCY5_NcadRFP_CPFundiff_DAPI_RFP_CY5_CFP\maxProjections';
ff = dir(maxproj_dir);
all_h5 = struct;
membrane_mask = struct;
mask_out_h5 = struct;
mask_out = struct;
q = 1;
for ii =1:size(ff,1)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,'_CY5_Probabilities.h5'))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5(q).file = ff(ii).name; 
membrane_mask(q).raw =all_h5(q).file;
mask_out_h5(q).raw = h5read(membrane_mask(q).raw,'/exported_data');
membr_lbl = 2;
mask_out(q).mask = squeeze(mask_out_h5(q).raw(membr_lbl,:,:));
q = q+1;
end
end