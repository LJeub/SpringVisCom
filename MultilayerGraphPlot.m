function [h_nodes, h_edges, h_layers]=MultilayerGraphPlot(xy,A,varargin)
% [h_nodes,h_edges]=GRAPHPLOT(xy,A,options)
%
% This program was written to take in xy-coordinates for a multilayer network
% as well as an array of adjacency matrices for each layer and draws each 
% layer of the network with the same coordinates in a different plane in
% 3D-space. Coordinates can be computed using MultilayerSpringVisCom.m
% (optionally emphasisng community structure) or by using other functions
% on the aggregate network. This program accepts a number of key-value
% options to fine-tune the layout.
%
% Inputs:
%
% xy: matrix of xy-coordinates for each physical node in the multilayer network
%
% A: Cell array of adjacency matrices for each layer
%
% Outputs:
%
% h_nodes: vector of node handles
%
% h_edges: handle to edges (either a single patch object or individual
%	lines depending on colour options)
%
% h_layers: handle to layer plane patches
%
% Optional Variables: (can be provided as an options struct or key-value
% pairs)
%
% alpha: factor for edge strength
%
% scores: vector with weights for each node (e.g. communities), used to
%       assign node colors from the colormap given by nodecolors and nodecolorlim
%       If a matrix is given for scores, the column of the matrix are
%       interpreted as different aspects, and the code draws a pie-chart
%       for each node, showing a nodes share of each aspect. This can be
%       useful e.g. for visualising overlapping communities. This variant
%       ignores the shapes option.
%
% shapes: marker shapes for nodes, either a single shape, giving all nodes the
%     same shape, or a vector of shapes, one for each node. (possible shapes are
%     the same as for LineSpec)
%
% pointsize: size of nodes (defaults to 7) can be scalar, or a vector with
% one value for each node
%
% edgewidth: linewidth of edges
%
% nodecolors: colormap for coloring nodes (using the weights provided by
% scores) (defaults to distinguishable colors, with 0 drawn black)
%
% nodecolorlim: clim for nodecolors (default: [min(SCORES), max(SCORES)])
%
% edgecolors: colormap for coloring edges using edge weights (defaults to
% grayscale)
%
% edgecolorlim: clim for edgecolors (default: [min(edgeweight,0),
% max(edgeweight)])
%
% randedges: randomise edge weights: 0: actual egeweights, 1: uniformly random
%
% edgethreshold: ignore edges below this weight (after randomised
% edgeweights)
%
% axispadding: extra padding around coordinates
%
% layercolor: color for layer planes ([0.5,0.5,0.5])
%
% layeralpha: transparency level for layer planes (0.1)
%
% layerlabels: labels for layer planes ([])
%
% labelfont: font for layer labels ('Helvetica')
%
% labelrotation: rotation for layer labels (30)
%
% aspectratio: Aspect ratio to use to rescale coordinates ([1,1,1])
%
% view: View specification ([-15, 30])
%
% drawlayers: Bool option to switch of drawing layer planes (true)
%
% This code uses the distinguishable_colors function available from
% http://www.mathworks.co.uk/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors

% Version: 1.3
% Date: Wed  9 May 2018 15:25:32 CEST
% Author: Lucas Jeub
% Email: ljeub@iu.edu

layer_width=length(A{1});
n_layers=length(A);
N=layer_width*n_layers;

parseArgs=inputParser();
parseArgs.KeepUnmatched=true;
addParameter(parseArgs,'axispadding',0.1);
addParameter(parseArgs,'layercolor',[0.5,0.5,0.5]);
addParameter(parseArgs,'layeralpha',0.1);
addParameter(parseArgs,'layerlabels',[]);
addParameter(parseArgs,'labelfont','Helvetica');
addParameter(parseArgs,'labelrotation',30);
addParameter(parseArgs,'aspectratio',[1,1,1]);
addParameter(parseArgs,'view',[-15,30]);
addParameter(parseArgs,'scores',[]);
addParameter(parseArgs,'drawlayers',true);

parse(parseArgs,varargin{:});
options=parseArgs.Results;
options.isset=@(opt) ~isempty(options.(opt));
plot_options=parseArgs.Unmatched;

aspect_ratio=options.aspectratio;
options.axispadding=options.axispadding./aspect_ratio;

for i=1:n_layers
    xyl((1:layer_width)+(i-1)*layer_width,1)=i/aspect_ratio(1);
    xyl((1:layer_width)+(i-1)*layer_width,2)=xy(:,1)/aspect_ratio(2);
    xyl((1:layer_width)+(i-1)*layer_width,3)=xy(:,2)/aspect_ratio(3);
end

pos=0;
for i=1:n_layers
    [ii,jj,vv]=find(A{i});
    Ai((1:length(ii))+pos)=ii+(i-1)*layer_width;
    Aj((1:length(jj))+pos)=jj+(i-1)*layer_width;
    Av((1:length(vv))+pos)=vv;
    pos=pos+length(ii);
end

AS=sparse(Ai,Aj,Av,N,N);

k=sum(AS+AS');


%options.pointsize=options.pointsize(:);
%options.pointsize=options.pointsize*double(logical(k));

if iscell(options.scores)
    plot_options.scores=[];
    for i=1:length(options.scores)
        plot_options.scores=[plot_options.scores;options.scores{i}];
    end
else
    plot_options.scores=options.scores(:);
    if length(plot_options.scores)<N
        plot_options.scores=sparse(plot_options.scores,1,1,N,1);
    end
end

view(options.view);
xmin=1/aspect_ratio(1)-options.axispadding(1);
xmax=(n_layers)/aspect_ratio(1)+options.axispadding(1);
ymin=min(xyl(:,2));
ymin=(ymin-options.axispadding(2));
ymax=max(xyl(:,2));
ymax=(ymax+options.axispadding(2));
zmin=min(xyl(:,3));
zmin=(zmin-options.axispadding(3));
zmax=max(xyl(:,3));
zmax=(zmax+options.axispadding(3));
axis equal
axis([xmin,xmax,ymin,ymax,zmin,zmax]);
set(gca,'clipping','off')
hold on
[h_nodes_out, h_edges_out]=GraphPlot(xyl,AS,plot_options);

if options.drawlayers
    for i=1:n_layers
        patch_h(i)=patch([i,i,i,i]/aspect_ratio(1),[ymin,ymax,ymax,ymin],...
            [zmin,zmin,zmax,zmax],options.layercolor,'linestyle','none');
    end
    set(patch_h,'facealpha',options.layeralpha);
else
    patch_h=[];
end

axis equal
axis([xmin,xmax,ymin,ymax,zmin,zmax]);
set(gca,'clipping','off')
if options.isset('layerlabels')
    for i=1:length(options.layerlabels)
        text(i/aspect_ratio(1),ymax,zmin,options.layerlabels{i},....
            'HorizontalAlignment','left','VerticalAlignment','cap',...
            'rotation',options.labelrotation,'Fontname',options.labelfont);
    end
end

set(gca,'ydir','reverse')

if nargout>0
    h_nodes=h_nodes_out;
    h_edges=h_edges_out;
    h_layers=patch_h;
end

end

