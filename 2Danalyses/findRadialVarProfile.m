function [radVarMean,radVarStd] = findRadialVarProfile(gfpDat,rfpDat,binSz,varName)
%FINDRADIALVARPROFILE finds the average value of the given variable within each radial bin
%of the given microcolony dataset. Note that it is assumed that the positional data
%has already been centred (as should be performed by cleanPosData.m
%anyway). Data from the RFP and GFP channels are pooled.
%
%   INPUTS:
%       -gfpDat: Data for the GFP-labelled microcolonies, as output by
%       cleanPosData
%       -rfpDat: Data for the RFP-labelled microcolonies
%       -binSz: Desired bin width for the radial histogram
%       -varName: String indicating which of the fields of gfpDat/rfpDat
%       you want to generate the radial plot of.
%
%   OUTPUTS:
%       -radVarMean: The radial average of the specified variable (centred
%       on (0,0)).
%       -radVarStd: The radial standard deviation of the specified variable. 
%
%   Author: Oliver J. Meacock, (c) 2020

allCents = [gfpDat.pos;rfpDat.pos];
allVar = [gfpDat.(varName);rfpDat.(varName)];

%Bin and average variable based on radial distance from centre of print
radDists = sqrt(sum(allCents.^2,2));
[N,Edges,Bins] = histcounts(radDists,'BinWidth',binSz);

radVarMean = zeros(max(Bins),1);
radVarStd = zeros(max(Bins),1);
for i = 1:max(Bins)
    radVarMean(i) = mean(allVar(Bins==i));
    radVarStd(i) = std(allVar(Bins==i));
end