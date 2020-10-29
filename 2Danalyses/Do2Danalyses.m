clear all
close all

%This script allows automated analysis of microcolony print dimensions (x-
%and y- dimensions) and radial profiles of microcolony density and averaged
%nearest-neighbour distance.

root = 'D:\RavAna'; %Name of the location where all your data is stored
branches = {'7 x 8 x 1\s_vs_s_si080', '7 x 8 x 1\s_vs_s_si094', '8 x 8 x 1\s_vs_s_si067', '8 x 8 x 1\s_vs_s_si088'}; %List of input datasets (directories)
twigs = {'rfp', 'gfp'}; %Names of sub-directories containing different channel data for each dataset
leaf = 'CellFeatures.mat'; %Name of the file output by FAST containing the microcolony data (position and fluorescence intensity)

szsX = cell(size(branches));
szsY = cell(size(branches));

binSz = 25; %Spacing of bins in the radial direction
verbose = false;

for b = 1:size(branches,2)
    load([root,filesep,branches{b},filesep,twigs{1},filesep,leaf])
    featsRFP = trackableData;
    load([root,filesep,branches{b},filesep,twigs{2},filesep,leaf])
    featsGFP = trackableData;
    
    szX = zeros(size(featsRFP.Centroid));
    szY = zeros(size(featsGFP.Centroid));
    
    for p = 1:size(featsRFP.Centroid,1)        
        [gfpDat,rfpDat] = cleanPosData(featsRFP,featsGFP,p);
        
        allCents = [gfpDat.pos;rfpDat.pos];
        
        szX(p) = max(allCents(:,1))-min(allCents(:,1));
        szY(p) = max(allCents(:,2))-min(allCents(:,2));
        
        if b == 2 %Only run these analyses on 0.94 SI dataset
            radDensProfs{p} = findRadialDensityProfile(gfpDat,rfpDat,binSz);
            [radDistProfs{p},radDists{p},meanNNdists{p}] = findRadialDistanceProfile(gfpDat,rfpDat,binSz);
            [radArMeans{p},radArStds{p}] = findRadialVarProfile(gfpDat,rfpDat,binSz,'area');
            
            %Truncate radial profiles so you don't include bins that extend
            %outside the limits of the square (ish) print
            minDim = min(szX(p),szY(p))/2;
            maxBin = floor(minDim/binSz);
            radDensProfs{p} = radDensProfs{p}(1:maxBin);
            radDistProfs{p} = radDistProfs{p}(1:maxBin);
        end
        
        if verbose
            subplot(1,2,1)
            plot(featsRFP.Centroid{p}(:,1),featsRFP.Centroid{p}(:,2),'r.')
            hold on
            plot(featsGFP.Centroid{p}(:,1),featsGFP.Centroid{p}(:,2),'g.')
            hold off
            axis equal
            
            subplot(1,2,2)
            plot(rfpDat.pos(:,1),rfpDat.pos(:,2),'r.')
            hold on
            plot(gfpDat.pos(:,1),gfpDat.pos(:,2),'g.')
            hold off
            axis equal
            
            pause(1)
        end
    end
    
    szsX{b} = szX;
    szsY{b} = szY;
end

save('C:\Users\olijm\Desktop\RavAna\ColonyDimensions.mat','szsX','szsY')
save('C:\Users\olijm\Desktop\RavAna\RadialProfiles.mat','radDensProfs','radDistProfs','meanNNdists','radDists')