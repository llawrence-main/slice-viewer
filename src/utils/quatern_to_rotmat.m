function [R, R_qfac] = quatern_to_rotmat(quatern_bcd, qfac)
%{
Returns the rotation matrix R based on quaternion parameters b,c,d

IN
quatern_bcd (vector): 3-vector of quaternion parameters [b,c,d]
qfac (scalar): qfac parameter

OUT
R (matrix): rotation matrix (proper)
R_qfac (matrix): rotation matrix with qfac multiple through (potentially improper)

NOTES
reference: https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html
%}

% determine quaternion parameters
b = quatern_bcd(1); c=quatern_bcd(2); d=quatern_bcd(3);
a = sqrt(1.0-(b*b+c*c+d*d));

% calculate (proper) rotation matrix
r1 = [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ];
r2 = [ 2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ];
r3 = [ 2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ];
R = [r1; r2; r3];

% compute rotation matrix with qfac multiplied through (potentially
% improper)
qfac_mat = [...
    1,0,0;...
    0,1,0;...
    0,0,qfac];
R_qfac = R*qfac_mat;

end