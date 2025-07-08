%%% This script shows an example of loading a NIfTI image and a mask and
%%% displaying an axial slice of the image with the mask overlaid as a
%%% contour. The image is a T1-weighted MRI scan of a phantom comprising 
%%% three large bottles and nine small tubes, all filled with water, 
%%% submerged in a water bath. One of the large bottles is contoured in the
%%% mask.

%% add path
addpath(genpath('src'));

%% load data
filename_img = 'data/t1w_3mm.nii.gz';
filename_mask = 'data/t1w_3mm_mask.nii.gz';
nii = nii_tool('load', filename_img);
mask = nii_tool('img', filename_mask);

%% call slice viewer to view axial slice with overlaid contour
view_plane = 'axial';
slice_number = 80;
figure;
view_slice(nii, view_plane, slice_number,...
    'Contours', mask);