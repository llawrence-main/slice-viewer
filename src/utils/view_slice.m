function image_handle = view_slice(data,varargin)
%{
Displays an image slice and overlaid contours. See the comments next to the
parameter declarations below for a description of the additional arguments.

IN
data (multiple types): data to view: a NIfTI structure, a 3D volume, or a matrix

OUT
image_handle: handle to image

NOTES
- The typical call is "view_slice(data, view, slice)", where "view" is one of
{'sagittal', 'axial', 'coronal'} and "slice" is the slice number to show
%}

% Input parser -- positional arguments
iparser = inputParser;
view_validation =@(view) assert(any(strcmp(view,{'sagittal','coronal','axial'})),'view must be one of {sagittal,coronal,axial}');
addOptional(iparser,'View','',view_validation); % View: sagittal, coronal, or axial
addOptional(iparser,'SliceNumber',[]); % The slice number to show

% Input parser -- parameter arguments
addParameter(iparser,'VoxelDimensions',[]); % Dimensions of voxels to apply appropriate scaling to image (not needed if nifti structure passed as data)
addParameter(iparser,'ColorMap',[]); % Color map for imagesc
addParameter(iparser,'ColorLimits',[]); % Limits of color map
addParameter(iparser,'ShowLabels',false); % Show anatomical labels?
addParameter(iparser,'Contours',[]); % ROIs to show as contours as a 3D matrix (2D contours) or 4D matrix (3D contours)
addParameter(iparser,'ContourColors',[]); % Colors to use for contours
addParameter(iparser,'ContourStyles','-'); % Line style for contours
addParameter(iparser,'ContourLineWidths',2); % Line width for contours
addParameter(iparser,'ContourDimensions',3); % Dimensions of the ROIs in the stack of contours (2D or 3D)
addParameter(iparser,'ColorbarLabel',''); % Label for color bar
addParameter(iparser,'ContourType','curve'); % Type of contour display (curve or wash)
addParameter(iparser,'WashAlpha',0.3); % Transparency for color wash
addParameter(iparser,'GaussSmooth', []); % Apply Gaussian smoothing with given sigma in mm before display
addParameter(iparser,'IgnoreOblique',false); % Round rotation matrix to nearest integer (even if oblique volume)?
addParameter(iparser,'Overlay',[]); % structure containing fields {img=image to overlay, limits=2-vector containing lower and upper intensity limits, colormap=colormap to use for overlay, and clims=color limits to ues for overlay}
parse(iparser,varargin{:});

% Determine if data is a NIfTI struct, a volume, or a matrix
if isstruct(data) && isfield(data,'img') && isfield(data,'hdr')
    datatype = 'nii';
elseif ndims(data)==3
    datatype = 'vol';
elseif ismatrix(data)
    datatype = 'matrix';
else
    assert(false,'data must be a NIfTI structure, a volume, or a matrix.');
end

% Check inputs depending on type of data
if strcmp(datatype,'matrix')
    view = 'axial';
else
    view = iparser.Results.View;
    view_validation(view);
    slice_number = iparser.Results.SliceNumber;
    assert(~isempty(slice_number),'If data is a volume or NIfTI structure, then a slice number must be passed.');
end

% Set orientation of axes
ignore_oblique = iparser.Results.IgnoreOblique;
if strcmp(datatype,'nii')
    if isfield(data.hdr, 'quatern_bcd')
        quatern_bcd = data.hdr.quatern_bcd;
    else
        quatern_bcd = [data.hdr.quatern_b, data.hdr.quatern_c, data.hdr.quatern_d];
    end
    qfac = data.hdr.pixdim(1);
    [~, R_qfac] = quatern_to_rotmat(quatern_bcd, qfac);
    if ignore_oblique
        R_qfac = round(R_qfac);
    end
    [nz_dims, is_oblique] = findCardinalDims(R_qfac);
    if is_oblique
        warning('Rotation matrix is not diagonal; cannot determine orientation of axes. Display geometry may be incorrect.');
        ax_or = 'RAS';
    else
        ax_or_list = {...
            'L', 'R';...
            'P', 'A';...
            'I', 'S'};        
        ax_or = '';
        for ix_dim = 1:3
            ix_dir = R_qfac(ix_dim, nz_dims(ix_dim))/2 + 3/2; % map [-1, 1] to [1, 2]
            ax_or = strcat(ax_or,ax_or_list{ix_dim, ix_dir});
        end
    end
else
    ax_or = 'RAS';
end

% Define slice(s) to show
if strcmp(datatype,'matrix')
    slice = data;
else
    if strcmp(datatype,'nii')
        vol = double(data.img);
    else
        vol = data;
    end
    slice = subset_to_slice(vol, view, slice_number);   
end

% Define voxel dimensions
if strcmp(datatype,'nii')
    vox_dims = data.hdr.pixdim(2:4);
else
    vox_dims = iparser.Results.VoxelDimensions;
    if isempty(vox_dims)
        vox_dims = ones(1,3);
    end
end

% apply smoothing if requested
smooth_sigma_mm = iparser.Results.GaussSmooth;
if ~isempty(smooth_sigma_mm)
    % make smoothing sigma the correct size
    if length(smooth_sigma_mm)==1
        smooth_sigma_mm = smooth_sigma_mm*ones(1,2);
    end

    % extract appropriate voxel dimensions
    switch view
        case 'sagittal'
            vox_dims_inplane = vox_dims(2:3);
        case 'coronal'
            vox_dims_inplane = vox_dims([1,3]);
        otherwise
            vox_dims_inplane = vox_dims(1:2);
    end

    % apply smoothing
    smooth_sigma_vox = smooth_sigma_mm./vox_dims_inplane;
    slice = imgaussfilt(slice, smooth_sigma_vox);
end

% Plot slice
image_handle = imagesc(slice);

% Color map
cmap = iparser.Results.ColorMap;
if isempty(cmap)
    colormap(gca,gray(256));
else
    colormap(gca,cmap);    
end

% Color limits
clims = iparser.Results.ColorLimits;
if ~isempty(clims)
    caxis(clims);
end

% Colorbar label
cbar_label = iparser.Results.ColorbarLabel;
if ~isempty(cbar_label)
    hcbar = colorbar;
    ylabel(hcbar,cbar_label);
end

% Declare the axes
ax = gca;

% Add overlay
overlay = iparser.Results.Overlay;
if ~isempty(overlay)

    % mask only where the image lies between the limits
    mask = (overlay.limits(1) < overlay.img) & (overlay.img < overlay.limits(2));
    overlay.img(~mask(:)) = 0;

    % subset overlay to appropriate slice
    overlay.img = subset_to_slice(overlay.img, view, slice_number);
    mask = subset_to_slice(mask, view, slice_number);

    % create another axis and add overlay
    ax_overlay = axes;
    imagesc(ax_overlay, overlay.img, 'alphaData', mask);
    colormap(ax_overlay, overlay.colormap);
    caxis(ax_overlay, overlay.clims);

    % update axis properties
    ax_overlay.Visible = 'off';
    linkprop([ax, ax_overlay], {'Position', 'XDir', 'YDir', 'DataAspectRatio', 'XTick', 'YTick'});

end

% Contours
contours_passed = iparser.Results.Contours;
contour_dimensions = iparser.Results.ContourDimensions;
contour_styles = iparser.Results.ContourStyles;
contour_line_widths = iparser.Results.ContourLineWidths;
assert(any(contour_dimensions==[2,3]),'Contour dimensions must be 2 or 3.');

% axis for contours
if isempty(overlay)
    ax_contour = ax;
else
    ax_contour = ax_overlay;
end

% Display contours if nonempty
if ~isempty(contours_passed)    
    if isstruct(contours_passed)&&~(isfield(contours_passed,'img')||isfield(contours_passed,'hdr'))
        % If contours is a structure, with each field a contour
        contour_names = fieldnames(contours_passed);
        num_contours = length(contour_names);
        size_contours = size(contours_passed.(contour_names{1}));
        contours = false([size_contours,num_contours]);
        for contour_no = 1:num_contours
            contour_name = contour_names{contour_no};
            contour = contours_passed.(contour_name);
            if contour_dimensions == 3
                contours(:,:,:,contour_no) = contour;
            else
                contours(:,:,contour_no) = contour;
            end
        end        
    elseif isstruct(contours_passed)&&isfield(contours_passed,'img')&&isfield(contours_passed,'hdr')
        % if contours passed is a nifti
        contours = logical(contours_passed.img);
    else
        % if contours passed is a volume
        contours = contours_passed;        
    end
    if contour_dimensions == 3
        num_contours = size(contours,4);
        % Subset to slice
        contours_slice = subset_to_slice(contours, view, slice_number);
    else
        num_contours = size(contours,3);
        contours_slice = contours;
    end
    % Define contour colours
    contour_colors = iparser.Results.ContourColors;
    if isempty(contour_colors)
        contour_colors = prism(num_contours);
    end
    % Define contour styles
    if length(contour_styles)==1
        contour_styles_deal = contour_styles;
        contour_styles = cell(1,num_contours);
        [contour_styles{1:num_contours}] = deal(contour_styles_deal);
    end
    % Define contour line widths
    if length(contour_line_widths)==1
        contour_line_widths = contour_line_widths*ones(1,num_contours);
    end
    % Define contour types
    contour_types = iparser.Results.ContourType;
    if ischar(contour_types)
        contour_types = repmat({contour_types},1,num_contours);
    end
    check_types =@(x) any(strcmp(x,{'curve','wash'}));
    assert(all(cellfun(check_types,contour_types)),'ContourType must be one of {curve,wash}');
    if any(strcmp(contour_types,'wash'))
        alph = iparser.Results.WashAlpha;
    end
    % Draw contours
    hold on
    for contour_no = 1:num_contours
        contour = squeeze(contours_slice(:,:,contour_no));
        switch contour_types{contour_no}
            case 'curve'
                contour_boundary_array = bwboundaries(contour); % Trace boundaries with 4-connectivity
                num_blobs = length(contour_boundary_array);
                for blob_no = 1:num_blobs
                    contour_boundary = contour_boundary_array{blob_no};
                    plot(ax_contour, contour_boundary(:,2),contour_boundary(:,1),...
                        'Color',contour_colors(contour_no,:),...
                        'LineStyle',contour_styles{contour_no},...
                        'LineWidth',contour_line_widths(contour_no));
                end
                
            case 'wash'
                wash = repmat(permute(contour_colors(contour_no,:),[3,1,2]),size(contour,1),size(contour,2)); % repeat colour vector
                h = imshow(wash, 'Parent', ax_contour);
                alpha_data = double(contour)*alph;
                set(h,'AlphaData',alpha_data);
        end
                
    end
    hold off
end

% Axis orientation
if strcmp(view,'sagittal')    
    if ax_or(2)=='A'
        xdir_or = 'normal';
    elseif ax_or(2)=='P'
        xdir_or = 'reverse';
    end
    if ax_or(3)=='S'
        ydir_or = 'normal';
    elseif ax_or(3)=='I'
        ydir_or = 'reverse';
    end
    ax.XDir = xdir_or;
    ax.YDir = ydir_or;    
    anatomical_labels_x = {'P','A'};
    anatomical_labels_y = {'I','S'};
    daspect_vec = [vox_dims(3),vox_dims(2),1];
elseif strcmp(view,'coronal')
    if ax_or(1)=='R'
        xdir_or = 'reverse';
    elseif ax_or(1)=='L'
        xdir_or = 'normal';
    end
    if ax_or(3)=='S'
        ydir_or = 'normal';
    elseif ax_or(3)=='I'
        ydir_or = 'reverse';
    end
    ax.XDir = xdir_or;
    ax.YDir = ydir_or;
    anatomical_labels_x = {'R','L'};
    anatomical_labels_y = {'I','S'};
    daspect_vec = [vox_dims(3),vox_dims(1),1];
elseif strcmp(view,'axial')
    if ax_or(1)=='R'
        xdir_or = 'reverse';
    elseif ax_or(1)=='L'
        xdir_or = 'normal';
    end
    if ax_or(2)=='A'
        ydir_or = 'normal';
    elseif ax_or(2)=='P'
        ydir_or = 'reverse';
    end
    ax.XDir = xdir_or;
    ax.YDir = ydir_or;
    anatomical_labels_x = {'R','L'};
    anatomical_labels_y = {'P','A'};
    daspect_vec = [vox_dims(2),vox_dims(1),1];
end

% Axis scaling
daspect(daspect_vec);

% Anatomical labels
show_labels = iparser.Results.ShowLabels;
if show_labels
    xticks(ax, xlim);
    yticks(ax, ylim);
    xticklabels(ax, anatomical_labels_x);
    yticklabels(ax, anatomical_labels_y);
else
    xticks(ax, []);
    yticks(ax, []);
end
end

function [nz_dims, is_oblique] = findCardinalDims(R)
% determines the non-zero dimensions of a rotation matrix, or tells the
% caller that the matrix is oblique

nz_dims = NaN(1, 3);
is_oblique = false;
for ix = 1:3
    if nnz(R(ix, :))==1
        nz_dims(ix) = find(R(ix, :));
    else
        is_oblique = true;
    end
end
end

function data_slice = subset_to_slice(data, view, slice)
% subsets 3D or 4D data to the requested slice in appropriate view

view = validatestring(view, {'axial', 'coronal', 'sagittal'});
assert((ndims(data) == 3)|(ndims(data) == 4), 'data must be 3D or 4D');

switch view
    case 'sagittal'
        data_slice = squeeze(data(slice, :, :, :));
    case 'coronal'
        data_slice = squeeze(data(:, slice, :, :));
    otherwise
        data_slice = squeeze(data(:, :, slice, :));
end

data_slice = permute(data_slice, [2,1,3]);

end