function [pop1pos,pop2pos] = alignMicrocoloniesSingleBoundary(pop1pos,pop2pos)
%ALIGNMICROCOLONIESSINGLEBOUNDARY aligns and centres microcolony data such
%that population 1 of microcolonies are on the left-hand side of the image,
%while population 2 are on the right-hand side.
%
%Note that this function assumes a simple geometry of population 1 abutting
%population 2 at a single interface, like so: [1|2].
%   
%   INPUTS:
%       -pop1pos: The original positions of N microcolonies of the first
%       type, an Nx2 matrix with the first column denoting x-coordinates and
%       the second column denoting y-coordinates
%       -pop2pos: The original positions of microcolonies of the second
%       type
%
%   OUTPUTS:
%       -pop1pos, pop2pos: Updated versions of the input vectors following
%       print alignment.
%
%   Author: Oliver J. Meacock, (c) 2020

meanX = mean([pop1pos(:,1);pop2pos(:,1)],1);
meanY = mean([pop1pos(:,2);pop2pos(:,2)],1);
pop1pos = pop1pos - repmat([meanX,meanY],size(pop1pos,1),1);
pop2pos = pop2pos - repmat([meanX,meanY],size(pop2pos,1),1);

%Convert to polar coordinates
[~,GFPrho] = cart2pol(pop1pos(:,1),pop1pos(:,2));
[~,RFPrho] = cart2pol(pop2pos(:,1),pop2pos(:,2));

%Remove everything more than 600um away from (apparent) centre of print, then recalculate
%positions relative to new centre
maxRad = 600;
pop1pos(GFPrho > maxRad,:) = [];
pop2pos(RFPrho > maxRad,:) = [];

meanX = mean([pop1pos(:,1);pop2pos(:,1)],1);
meanY = mean([pop1pos(:,2);pop2pos(:,2)],1);
pop1pos = pop1pos - repmat([meanX,meanY],size(pop1pos,1),1);
pop2pos = pop2pos - repmat([meanX,meanY],size(pop2pos,1),1);

PCAableDat = pop1pos;

C = cov(PCAableDat);
[V,D] = eig(C);

rot = [V(1,1),-V(1,2);V(1,2),V(1,1)];

pop1pos = pop1pos*rot;
pop2pos = pop2pos*rot;

%Check to see if GFP spots are generally in negative x region - rotate
%everything by 180 degrees if not.
if ~(mean(pop1pos(:,1)) < 0)
    pop1pos = pop1pos*[-1,0;0,1];
    pop2pos = pop2pos*[-1,0;0,1];
end