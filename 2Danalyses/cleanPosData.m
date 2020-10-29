function [gfpData,rfpData] = cleanPosData(featsRFP,featsGFP,p)
%CLEANPOSDATA removes microcolonies far away from the main print, removes
%whole-print average position and rotates into a basis such that print is
%in a nice x-y coordinate space.
%
%   INPUTS:
%       -featsRFP: Output of FAST's feature detection module (the
%       cellFeatures.mat file) for the RFP-labelled microcolonies.
%       -featsGFP: Likewise for GFP-labelled microcolonies.
%       -p: Allows specification of a frame in featsRFP and featsGFP, if
%       many microcolony images have been analysed as a stack.
%
%   OUTPUTS:
%       -gfpData: Cleaned microcolony data for GFP-labelled microcolonies, 
%       with distant microcolonies removed and microcolony positions
%       registered to the quadrilateral shape of the print. Fields are pos
%       (Nx2 vector denoting microcolony coordinates), area (Nx1 vector
%       denoting microcolony area) and mean (Nx1 vector denoting mean
%       intensity).
%       -rfpData: Likewise for RFP-labelled microcolonies.
%
%   Author: Oliver J. Meacock, (c) 2020

%% Part 1: Remove outlying microcolonies

RFPpos = featsRFP.Centroid{p};
GFPpos = featsGFP.Centroid{p};

allPos = [GFPpos;RFPpos];

initCent = nanmean(allPos,1); %First go through (before removing outlying microcolonies)
allPos = allPos-repmat(initCent,size(allPos,1),1);

radDist = sqrt(sum(allPos.^2,2));

%For a square sampled with two uniform distributions, the average point is
%L/6(sqrt(2)+log(1+sqrt2))) units away from the origin (assuming the square
%is centred on the origin). L is the side length.

%Can use this fact to estimate L from the radial
%distance distribution and so guesstimate the maximal distance a
%microcolony should be from the print centre to be included in downstream
%analyses.

distMean = nanmean(radDist);
Lguess = 6*distMean/(sqrt(2)+log(1+sqrt(2)));

badDists = radDist > 0.75*Lguess; %Formally, maximum distance should be from centre to corner of squrare (i.e. sqrt(2)/2*L). But this keeps life simple, and accounts for the possibility of microcolony density being greater in the middle of the colony.
badDistsGFP = badDists(1:size(GFPpos,1));
badDistsRFP = badDists(size(GFPpos,1)+1:end);

badDistsGFP = or(badDistsGFP,isnan(GFPpos(:,1)));%Further - get rid of any nan values that pop up due to microcolonies being too small
badDistsRFP = or(badDistsRFP,isnan(RFPpos(:,1)));

%% Part 2: Perform the actual transformations based on the cleaned microcolony data
cleanRFPpos = RFPpos(~badDistsRFP,:);
cleanGFPpos = GFPpos(~badDistsGFP,:);

cleanAllPos = [cleanRFPpos;cleanGFPpos];

newCent = nanmean(cleanAllPos,1);
cleanAllPos = cleanAllPos-repmat(newCent,size(cleanAllPos,1),1);
cleanRFPpos = cleanRFPpos-repmat(newCent,size(cleanRFPpos,1),1);
cleanGFPpos = cleanGFPpos-repmat(newCent,size(cleanGFPpos,1),1);

%This is a method inspired by the hexatic order parameter/ Fourier
%transforms
outerPts = convhull(cleanAllPos(:,1),cleanAllPos(:,2));
[thetOut,rhoOut] = cart2pol(cleanAllPos(outerPts,1),cleanAllPos(outerPts,2));
rhoQ = interp1(thetOut(2:end),rhoOut(2:end),linspace(-pi,pi,100),'linear','extrap'); %Should ensure number of points in resampling is divisible by 4.
rhoQ = rhoQ - mean(rhoQ);
cplxSum = sum(rhoQ.*exp(((1i*2*pi)/size(rhoQ,2))*(0:(size(rhoQ,2)-1))*4)); % This is essentially a DFT for k = 4 (the number of corners in a square)
rotAng = (pi-angle(cplxSum))/4; %Using the phase information in the Fourier representation to estimate the orientation of the square.

rot = [cos(rotAng),-sin(rotAng);sin(rotAng),cos(rotAng)];

gfpData.pos = cleanGFPpos*rot';
rfpData.pos = cleanRFPpos*rot';

gfpData.area = featsGFP.Area{p}(~badDistsGFP);
gfpData.mean = featsGFP.ChannelMean{p}(~badDistsGFP);

rfpData.area = featsRFP.Area{p}(~badDistsRFP);
rfpData.mean = featsRFP.ChannelMean{p}(~badDistsRFP);