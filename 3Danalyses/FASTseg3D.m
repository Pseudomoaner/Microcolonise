function outSeg = FASTseg3D(im,segSets)
%FASTSEG3D runs a similar segmentation routine to FAST, but generalised for
%three-dimensional imaging data. For details on the 'ridge detection'
%routine used, see:
%   -Lindeberg (1998): Edge Detection and Ridge Detection with Automatic
%   Scale Selection.
%   -Jin et al. (2013): Vascular tree segmentation in medical images using 
%   Hessian-based multiscale filtering and level set method.
%
%   INPUTS:
%       -im: The original 3-dimensional matrix denoting the input image.
%       -segSets: Structure defining the settings used in the segmentation.
%       Should contain the following fields:
%           -ridgeThresh: Threshold applied to the largest normalized
%           eigenvalue of the Hessian at each voxel; voxels with values
%           above this value will be marked true.
%           -ridgeScale: Width of the Gaussian smoothing kernal
%           -bumpThresh: The intensity-based binarization threshold
%           -ridgeErosion: Size of the structure element used to open the
%           binarized 3D image.
%           -minVol: Minimum volume of a microcolony to be detected (in
%           voxels)
%
%   OUTPUTS:
%       -outSeg: Three dimensional matrix of the same dimensions as im,
%       with detected objects marked with true voxels and background marked
%       as false voxels.
%
%   Author: Oliver J. Meacock, (c) 2020

%Set parameters
MedFiltSize = 3;

%Apply a median filter
im = medfilt3(im,[MedFiltSize,MedFiltSize,MedFiltSize]);

%Do feature scaling
im = (im - min(im(:)))/(max(im(:)) - min(im(:)));

%Threshold directly
Bumps = imbinarize(im,segSets.bumpThresh);

%Do ridge-detection segmentation
Ridges = bwRidgeCenterMod(im,segSets.ridgeScale,segSets.ridgeThresh);
se = strel('sphere',segSets.ridgeErosion);
Ridges = imerode(Ridges,se);
Ridges = imdilate(Ridges,se);

tempImg = and(Bumps, ~Ridges);

%Remove any holes that might have got into the image at this point
fillImg = imfill(tempImg,'holes');

%Do an image opening to clean microcolony boundaries and get rid of 'dust'
fillImg = imerode(fillImg,se);
fillImg = imdilate(fillImg,se);

fillImg = bwareaopen(fillImg,segSets.minVol);

outSeg = fillImg;