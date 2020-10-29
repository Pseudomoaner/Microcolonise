clear all
close all

%This script runs three dimensional segmentation by combining a Hessian
%eigenvalue measurement approach with a simple intensity threshold. 
%For more details on the Hessian eigenvalue approach, see Lindeberg (1998): 
%Edge Detection and Ridge Detection with Automatic Scale Selection (among
%others)

root = 'C:\Users\olijm\Desktop\RavAna\3D'; %Location of directory where all the data is located
branch = 'Colony'; %Generic name of the data subdirectories (will have an integer added onto the end in the main loop - i.e. becoming Colony1, Colony2 etc.)
RFPleaf = 'RFP.tif'; %Names of the two channel images in each subdirectory
GFPleaf = 'GFP.tif';

segLeaf = 'Segmentations.mat';
outLeaf = 'pixelCounts.mat';

noBranches = 6;

%3D segmentation parameters
segSettings.ridgeThresh = 0.05;
segSettings.ridgeScale = 15; %Width of the Gaussian smoothing kernal
segSettings.bumpThresh = 0.3; %The intensity-based binarization threshold
segSettings.ridgeErosion = 1;
segSettings.minVol = 30; %Minimum volume of a microcolony to be detected

GFPpx = zeros(noBranches,1);
RFPpx = zeros(noBranches,1);

for b = 1:noBranches
    rfpPath = [root,filesep,branch,num2str(b),filesep,RFPleaf];
    gfpPath = [root,filesep,branch,num2str(b),filesep,GFPleaf];
    
    %Get image dimensions
    info = imfinfo(rfpPath); %I'll assume the RFP and GFP images are the same dimensions
    noX = info(1).Width;
    noY = info(1).Height;
    noZ = size(info,1);
    
    fullGFP = zeros(noY,noX,noZ);
    fullRFP = zeros(noY,noX,noZ);
    for i = 1:noZ
        fullRFP(:,:,i) = imread(rfpPath,i);
        fullGFP(:,:,i) = imread(gfpPath,i);
    end
    
    segRFP = FASTseg3D(fullRFP,segSettings);
    segGFP = FASTseg3D(fullGFP,segSettings);
    
    RFPpx(b) = sum(segRFP(:));
    GFPpx(b) = sum(segGFP(:));

    save([root,filesep,branch,num2str(b),filesep,segLeaf],'segRFP','segGFP'); %Save the full 3D segmentations
end

save([root,filesep,outLeaf],'RFPpx','GFPpx') %Save the summary statistics