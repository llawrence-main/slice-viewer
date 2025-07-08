function slice_tgt = same_slice(nii_src,slice_src,nii_tgt,varargin)
%{
Returns the closest corresponding slice to the source slice in the target

IN
nii_src (struct): source nifti
slice_src (int): slice in source nifti
nii_tgt (struct): target nifti
plane (optional, str): plane of slice (axial, coronal, sagittal; default = axial)

OUT
slice_tgt (int): slice in target nifti
%}

% parse arguments
if nargin == 3
    plane = 'axial';
else
    plane = varargin{1};
end
assert(any(strcmp(plane,{'axial','coronal','sagittal'})),'plane is not valid');

% declare parameters for given plane
switch plane
    case 'axial'
        f = 'qoffset_z';
        idx = 4;
    case 'coronal'
        f = 'qoffset_y';
        idx = 3;
    case 'sagittal'
        f = 'qoffset_x';
        idx = 2;
end
        
% compute corresponding slice
slice_src = slice_src - 1; % switch from MATLAB 1-indexing to NIfTI 0-indexing
ax_offset_src = nii_src.hdr.(f);
ax_offset_tgt = nii_tgt.hdr.(f);
delta_ax_src = nii_src.hdr.pixdim(idx);
delta_ax_tgt = nii_tgt.hdr.pixdim(idx);
slice_tgt = (ax_offset_src-ax_offset_tgt+slice_src*delta_ax_src)/delta_ax_tgt;
slice_tgt = round(slice_tgt);
slice_tgt = slice_tgt + 1; % switch from NIfTI 0-indexing to MATLAB 1-indexing

end