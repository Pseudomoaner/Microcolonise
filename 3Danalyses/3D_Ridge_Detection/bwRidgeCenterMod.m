function bwridge = bwRidgeCenterMod(I,  scales, valuethresh)

%% extract the centerline of scale-space valleys
M_mag = im_scalablehess3( I, scales, 'dark');

if isempty(valuethresh)
    valuethresh = graythresh( M_mag );
end
bwridge = imbinarize( M_mag, valuethresh );
