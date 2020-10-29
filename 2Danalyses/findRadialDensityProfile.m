function radDensProf = findRadialDensityProfile(gfpDat,rfpDat,binSz)
%FINDRADIALDENSITYPROFILE finds the density of microcolonies within each radial bin
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
%       -radDensProf: The radially averaged density profile for the pooled 
%       GFP and RFP microcolony data.
%
%   Author: Oliver J. Meacock, (c) 2020

allCents = [gfpDat.pos;rfpDat.pos];

radDists = sqrt(sum(allCents.^2,2));
[N,Edges] = histcounts(radDists,'BinWidth',binSz);

%Find the area of each annular section
annArs = zeros(size(Edges,2)-1,1);
for i = 1:size(Edges,2)-1
    innerCirc = pi*(Edges(i).^2);
    outerCirc = pi*(Edges(i+1).^2);
    
    annArs(i) = outerCirc-innerCirc;
end

radDensProf = N'./annArs;