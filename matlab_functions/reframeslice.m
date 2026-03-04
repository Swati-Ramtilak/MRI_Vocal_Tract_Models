function imgout = reframeslice(img, nsl_row, nsl_col, st_slice, plane, step)
%function imgout  = reframeslice(inmtx, nsl_row, nsl_col, st_slice)
% inmtx           : input 3D matrix
% nsl_row         : # of slice (whole image) per column
% nsl_col         : # of slice (whole image) per row
% st_slice        : starting slice to be reframed
%plane            : 0 (axial), 1(coronal), 2(sagittal)
% step            : step size in slice direction
% imgout          : reframed image

if nargin < 6
    step = 1;
end
if nargin < 5
    plane = 0;
end
[nslices, nphases, imsize] = size(img);
imgout = zeros(nsl_row*nphases, nsl_col*imsize);
st_slice = max(0, st_slice);
end_slice = min(nsl_col*nsl_row*step+st_slice-1, nslices);
for k = st_slice:step:end_slice
    for m = 1:nphases
        row = floor((k-st_slice)/nsl_col/step) * nphases + m;
        for n = 1:imsize
            col = mod((k-st_slice)/step, nsl_col) * imsize + n;
            if plane == 0
                imgout(row,col) = img(k, m, n);
            elseif plane == 1
                imgout(row,col) = img(m, k, n);
            else
                imgout(row,col) = img(m, n, k);
            end
        end
    end
end
end