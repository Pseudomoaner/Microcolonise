function [radDistProf,radDists,meanDists] = findRadialDistanceProfile(gfpDat,rfpDat,binSz)
%FINDRADIALDISTANCEPROFILE finds the drop-drop spacing of microcolonies within each radial bin
%of the given microcolony dataset. Note that it is assumed that the positional data
%has already been centred (as should be performed by cleanPosData.m
%anyway).
%
%   INPUTS:
%       -gfpDat: Data for the GFP-labelled microcolonies, as output by
%       cleanPosData
%       -rfpDat: Data for the RFP-labelled microcolonies
%       -binSz: Desired bin width for the radial histogram
%
%   OUTPUTS:
%       -radDistProf: The radial average of the average distance beween
%       each microcolony and its nearest neighbours (defined with Delaunay
%       triangulation)
%       -radDists: The distance of each microcolony and the print center
%       (position (0,0))
%       -meanDists: The average distance beween each microcolony and its 
%       nearest neighbours (defined with Delaunay triangulation) 
%
%   Author: Oliver J. Meacock, (c) 2020

allCents = [gfpDat.pos;rfpDat.pos];

%% Part 1: Create a list of average distances between each microcolony and its (Delaunay triangulated) neighbours.
tris = delaunay(allCents(:,1),allCents(:,2));
meanDists = zeros(size(allCents,1),1);
for i = 1:size(allCents,1) %For each microcolony
    currTris = logical(sum(tris == i,2));
    currPts = unique(tris(currTris,:));
    currPts(currPts == i) = []; %Remove self from list
    dists = sqrt(sum((allCents(currPts,:)-repmat(allCents(i,:),size(currPts,1),1)).^2,2));
    meanDists(i) = mean(dists);
end

%% Part 2: Bin and average list based on radial distance from centre of print
radDists = sqrt(sum(allCents.^2,2));
[N,Edges,Bins] = histcounts(radDists,'BinWidth',binSz);

radDistProf = zeros(max(Bins),1);
for i = 1:max(Bins)
    radDistProf(i) = mean(meanDists(Bins==i));
end