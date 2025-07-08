function roi_slices = nonzero_slices(roi)
%{
Returns the axial slices over which an ROI is nonzero

IN
roi (3D volume): ROI array

OUT
slices (vector): slices over which ROI is not zero
%}

roi_flat = permute(roi,[3,1,2]);
roi_slices = find(sum(roi_flat(:,:),2));

end