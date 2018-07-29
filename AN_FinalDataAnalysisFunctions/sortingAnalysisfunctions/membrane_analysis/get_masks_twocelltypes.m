function [mask1_out_h5,mask2_out_h5,mask_out1,mask_out2,all_h5_type1,all_h5_type2]=get_masks_twocelltypes(maxproj_dir,str_type1,str_type2)

ff = dir(maxproj_dir);
all_h5_type1 = struct;
all_h5_type2 = struct;
nuc_mask1 = struct;% for cell type 1
nuc_mask2 = struct;% for cell type 2
mask1_out_h5 = struct;
mask_out1 = struct;
mask2_out_h5 = struct;
mask_out2 = struct;
q = 1;
s = 1;
for ii =1:size(ff,1)    
    if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,str_type1))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5_type1(q).file = ff(ii).name; 
nuc_mask1(q).raw =all_h5_type1(q).file;
mask1_out_h5(q).raw = h5read(nuc_mask1(q).raw,'/exported_data');
nuc_lbl = 2;
mask_out1(q).mask = squeeze(mask1_out_h5(q).raw(nuc_lbl,:,:));
q = q+1;
    end
if ~isdir(ff(ii).name) && ~isempty(strfind(ff(ii).name,str_type2))
        disp(['loading: ' num2str(ff(ii).name)])
        all_h5_type2(s).file = ff(ii).name; 
nuc_mask2(s).raw =all_h5_type2(s).file;
mask2_out_h5(s).raw = h5read(nuc_mask2(s).raw,'/exported_data');
nuc_lbl = 2;
mask_out2(s).mask = squeeze(mask2_out_h5(s).raw(nuc_lbl,:,:));
s = s+1;
    end
end

end