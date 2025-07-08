function [xlims,ylims] = zoom_limits(roi,margin,varargin)
%{
Determines the x and y axis limits required to zoom in on largest
connected component of an ROI with a given margin

IN
roi (3D volume): ROI array
margin (int): margin in number of voxels around the ROI bounding box
view (char, optional): anatomical view (axial, coronal, sagittal)
(default: axial)
slice (int, optional): slice to restrict to (default: no restriction)

OUT
xlims (2-vector): x-axis limits
ylims (2-vector): y-axis limits
%}

% parse arguments
if nargin>3
    view = varargin{1};
    slice = varargin{2};
elseif nargin>2
    view = varargin{1};
    slice = [];
else
    view = 'axial';
    slice = [];
end

if ~isempty(slice)
    roi_tmp = false(size(roi));
    switch view
        case 'sagittal'
            roi_tmp(slice,:,:) = roi(slice,:,:);
        case 'coronal'
            roi_tmp(:,slice,:) = roi(:,slice,:);
        case 'axial'            
            roi_tmp(:,:,slice) = roi(:,:,slice);        
    end
    roi = roi_tmp;
end

% determine largest connected component
cc = bwconncomp(roi);
[~,idx] = max(cellfun(@numel,cc.PixelIdxList));

% compute region properties
rp = regionprops(cc);

% extract bounding box of largest connected component
bb = rp(idx).BoundingBox;

% compute i,j,k-axis limits
imin = ceil(bb(2)) - margin;
imax = imin + bb(5) + 2*margin - 1;
jmin = ceil(bb(1)) - margin;
jmax = jmin + bb(4) + 2*margin - 1;
kmin = ceil(bb(3)) - margin;
kmax = kmin + bb(6) + 2*margin - 1;

% clip to stay within frame of image
sz = size(roi);
imin = max(imin, 1); jmin = max(jmin, 1); kmin = max(kmin, 1);
imax = min(imax, sz(1)); jmax = min(jmax, sz(2)); kmax = min(kmax, sz(3));

% get x and y axis limits based on view
switch view
    case 'sagittal'
        xlims = [jmin,jmax];
        ylims = [kmin,kmax];
    case 'coronal'
        xlims = [imin,imax];
        ylims = [kmin,kmax];
    case 'axial'
        xlims = [imin,imax];
        ylims = [jmin,jmax];
end

end