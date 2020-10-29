clear all
close all

%This script creates plots of the radial distributions of fluorescently
%labeled susceptible and resistant strains when grown next to each other.
%Also plots the result of a susceptible/susceptible control.

root = 'K:\RavAna';
branches = {'s_vs_p', 's_vs_s_2'};
twigs = {'rfpchannel', 'gfpchannel'};
leaf = 'CellFeatures.mat';

pop1Tags = {'P_{sp}','S_{ss}'};
pop2Tags = {'S_{sp}','S_{ss}'};

binEdges = 25:50:475;

GFPMeans = cell(size(branches));
RFPMeans = cell(size(branches));

for br = 1:size(branches,2)
    load([root,filesep,branches{br},filesep,twigs{1},filesep,leaf],'trackableData')
    RFPdat = trackableData;
    load([root,filesep,branches{br},filesep,twigs{2},filesep,leaf],'trackableData')
    GFPdat = trackableData;
    
    for i = 1:size(GFPdat.Centroid,1)        
        GFPpos = GFPdat.Centroid{i};
        RFPpos = RFPdat.Centroid{i};
        GFParea = GFPdat.Area{i};
        RFParea = RFPdat.Area{i};
        
        GFPpos(isnan(GFParea),:) = [];
        RFPpos(isnan(RFParea),:) = [];
        GFParea(isnan(GFParea),:) = [];
        RFParea(isnan(RFParea),:) = [];
        
        [GFPpos,RFPpos] = alignMicrocoloniesSingleBoundary(GFPpos,RFPpos);
        
        %Reconvert to polar coordinates
        [GFPthet,GFPrho] = cart2pol(GFPpos(:,1),GFPpos(:,2));
        [RFPthet,RFPrho] = cart2pol(RFPpos(:,1),RFPpos(:,2));
        
        %Bin by radii and get means
        [~,~,GFPbin] = histcounts(GFPrho,binEdges);
        [~,~,RFPbin] = histcounts(RFPrho,binEdges);
        
        for j = 1:size(binEdges,2)-1
            GFPMean(i,j) = mean(GFParea(GFPbin == j));
            RFPMean(i,j) = mean(RFParea(RFPbin == j));
        end
    end
    GFPMeans{br} = GFPMean;
    RFPMeans{br} = RFPMean;
end

fh = figure(1);
hold on

options.x_axis = binEdges(1:end-1) + diff(binEdges)/2;
options.error = 'std';
options.alpha = 0.3;
options.line_width = 2;
options.fHandle = fh;

for i = 1:size(branches,2)
    for j = 1:size(branches,2)
        options.axHandle = subplot(size(branches,2),size(branches,2),(i-1)*size(branches,2)+j);
        options.color_area = [0,1,0];
        options.color_line = [0,1,0];
        
        hsG = plot_areaerrorbar(GFPMeans{i},options);
        
        options.color_area = [1,0,0];
        options.color_line = [1,0,0];
        
        hsR = plot_areaerrorbar(RFPMeans{j},options);
        
        xlabel('Radial distance (um)')
        ylabel('Microcolony area (um^2)')
        
        options.axHandle.LineWidth = 1;
        
        legend([hsR(2),hsG(2)],{pop1Tags{j},pop2Tags{i}},'Location','SouthEast')
        axis([50,425,0,250])
    end
end

% hold on
% for i = 1:size(GFPMean,1)
%     plot(options.x_axis,GFPMean(i,:),'.','MarkerEdgeColor',[0,0.7,0])
%     plot(options.x_axis,RFPMean(i,:),'.','MarkerEdgeColor',[0.7,0,0])
% end