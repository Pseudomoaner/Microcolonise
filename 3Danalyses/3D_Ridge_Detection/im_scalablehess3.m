function Magnitude = im_scalablehess3(im, scale_range, ridge_property)

if (nargin == 2)
    ridge_property = 'dark';
end

dim = 0;
for scale = scale_range
    dim = dim + 1;
    Magnitude(:,:,:,dim) = im_hessangle3(im, scale);
end

if strcmp(ridge_property, 'bright')
    Magnitude = -Magnitude;
elseif strcmp(ridge_property, 'dark')
    ;
else
    error('ridge_property error @ scalablehess2');  
end

Magnitude(Magnitude<0) = 0;
Magnitude = (Magnitude-min(Magnitude(:)))/(max(Magnitude(:))-min(Magnitude(:)) + realmin); %Set to vary between 0 and 1
