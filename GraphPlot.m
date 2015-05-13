function [h_nodes_out,h_edges_out]=GraphPlot(xy,W,varargin)
% [h_nodes,h_edges]=GRAPHPLOT(xy,W,options)
%
% This program was written to take in xy-coordinates of a network as well as
% the adjacency matrix W, to graph the network using these coordinates.  We
% have many programs that can calculate these coordinates either respecting
% community structure (SpringVisCom.m, KamadaC.m, fruc_reinC.m, KKFR.m, FRKK.m) or ignoring
% it (Kamada.m, fruc_rein.m).  This program accepts a number of key-value
% options to fine-tune the layout.
%
% Inputs:
%
% xy: matrix of xy-coordinates for each node in the network
%
% W: Adjacency Matrix for the network
%
% Outputs:
%
% h_nodes: vector of node handles
% h_edges: handle to edges (either a single patch object or individual
%	lines depending on colour options)
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
%
% Alternatively use function call:
%
% [h_nodes,h_edges]=GRAPHPLOT2D(xy,W,alpha,scores,shapes,pointsize,...)
% where the order of options is as given above.
%
% This code uses the distinguishable_colors function available from
% http://www.mathworks.co.uk/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors

% Version: 1.2
% Date: Tue 13 May 2014 17:03:19 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk



%% Clean input arguments
% number of nodes
N=length(W);
base_radius=1;
AngularResolution=30;

% check inputs
if size(xy,1)~=N
    if size(xy,2)==N
        xy=xy';
    else
        size        error('coordinates are of the wrong size for provided adjacency matrix')
    end
end

% transform adjacency matrix into edge-list
[indrow,indcol,weight]=find(W);
if isequal(W,W') % undirected
    ind=find(indcol>indrow);
    edges=[indcol(ind),indrow(ind),weight(ind)];
    directed=false;
else
    edges=[indcol,indrow,weight];
    directed=true;
end

% check if hold is on
is_hold=ishold;

%% Set up options

% default options
options=OptionStruct('alpha',1,'scores',ones(N,1),'shapes','.','pointsize',7,'edgewidth',1,'nodecolors',[],'nodecolorlim',[],...
    'edgecolors',[],'edgecolorlim',[],'randedges',0,'edgethreshold',[],'legendlabels',[],'baseradius',base_radius);

% parse options
options.set(varargin);

base_radius=options.baseradius;
% if ~isempty(varargin)
%     if ischar(varargin{1})||isstruct(varargin{1})
%         options.set(varargin);
%     else
%         % assume old form of function given (backwards compatibility)
%         names=options.options;
%         input=struct([]);
%         for i=1:length(varargin)
%             input(1).(names{i})=varargin{i};
%         end
%         options.set(input);
%     end
% end
% make scores full matrix
options.scores=full(options.scores);
if isvector(options.scores)
    options.scores=options.scores(:);
end

if ~isempty(edges)
    % optionally randomise edgeweights
    edges(:,3)=(1-options.randedges)*edges(:,3)+max(edges(:,3))*rand(size(edges,1),1)*options.randedges;
    if options.isset('edgethreshold')
        edges=edges(edges(:,3)>options.edgethreshold,:);
    end
    
    
    %sort edges to plot smallest weight first
    [~,s]=sort(edges(:,3));
    edges=edges(s,:);
end
% parse named colors
if ischar(options.nodecolors)||iscellstr(options.nodecolors)
    options.nodecolors=colorstr2rgb(options.nodecolors);
end

if ischar(options.edgecolors)||iscellstr(options.edgecolors)
    options.edgecolors=colorstr2rgb(options.edgecolors);
end

% set up nodecolorlim
if ~isset(options,'nodecolorlim')
    if size(options.scores,2)==1
        options.nodecolorlim=[min(min(options.scores),0),max(options.scores)];
    else
        options.nodecolorlim=[0,size(options.scores,2)];
    end
end

% set up edgecolorlim
if ~isempty(edges)
    if ~isset(options,'edgecolorlim')
        options.edgecolorlim=[min(min(edges(:,3)),0),max(edges(:,3))];
    end
else
    options.edgecolorlim=[0,1];
end

% set up node shapes
if length(options.shapes)==1
    shapes=repmat(options.shapes,N,1);
else
    shapes=options.shapes(:);
end

% set up point sizes
if length(options.pointsize)==1
    point_size=ones(N,1)*options.pointsize;
else
    point_size=options.pointsize;
end

% set up edgewidth
edgewidth=options.edgewidth;


%% set up edge colors and edge plotting function
if edgewidth>0
    
    if isset(options,'edgecolors')
        edgemap=options.edgecolors;
        nedgemap=size(edgemap,1);
        if nedgemap>1
            if options.alpha~=1
                edgemap=interp1(linspace(0,1,nedgemap).^options.alpha,edgemap,linspace(0,1,nedgemap),'pchip');
            end
            u_weight=unique(edges(:,3));
            if length(u_weight)==size(edgemap,1)
                edgecolor = @(weight) edgemap(u_weight==weight,:);
            else
                edgecolor = @(weight) interp1(linspace(options.edgecolorlim(1),options.edgecolorlim(2),size(edgemap,1)),edgemap,weight,'nearest');
                
            end
        else
            edgecolor=@(weight) edgemap;
        end
        
        switch size(xy,2)
            case 2
                plot_edges = @plot_edges_fixed_color_2;
            case 3
                plot_edges = @plot_edges_fixed_color_3;
            otherwise
                error('need 2 or 3 dimensional coordinates');
        end
    else
        a=linspace(0.9,.2,size(colormap,1));
        edgemap=repmat(a',1,3);
        if options.alpha~=1
            edgemap=interp1(linspace(0,1,50).^options.alpha,edgemap,linspace(0,1,50),'pchip');
        end
        if directed
            edgecolor=@(weight) weight;
            switch size(xy,2)
                case 2
                    plot_edges =@plot_edges_fixed_color_2;
                case 3
                    plot_edges = @plot_edges_fixed_color_3;
                otherwise
                    error('need 2 or 3 dimensional coordinates');
            end
        else
            
            edgecolor=@(weight) repmat(weight(:)',2,1);
            switch size(xy,2)
                case 2
                    plot_edges=@() patch(reshape(xy(edges(:,1:2),1),size(edges,1),2)',reshape(xy(edges(:,1:2),2),size(edges,1),2)',edgecolor(edges(:,3)),'linewidth',edgewidth,'FaceColor','none','EdgeColor','flat');
                case 3
                    plot_edges=@() patch(reshape(xy(edges(:,1:2),1),size(edges,1),2)',reshape(xy(edges(:,1:2),2),size(edges,1),2)',reshape(xy(edges(:,1:2),3),size(edges,1),2)',edgecolor(edges(:,3)),'linewidth',edgewidth,'FaceColor','none','EdgeColor','flat');
                otherwise
                    error('need 2 or 3 dimensional coordinates');
            end
        end
    end
    
    if ~is_hold
        colormap(edgemap);
    end
else
    plot_edges=@() [];
end

    function h=plot_edges_fixed_color_2
        h=zeros(size(edges,1),1);
        if directed
            for e=1:size(edges,1)
                h(e)=arrow(xy(edges(e,1:2),1),xy(edges(e,1:2),2),[0,0],[0,0,1],edgecolor(edges(e,3)),edgewidth*base_radius,point_size(edges(e,2))*base_radius*1.01);
            end
        else
            for e=1:size(edges,1)
                h(e)=plot(xy(edges(e,1:2),1),xy(edges(e,1:2),2),'linewidth',edgewidth,'color',edgecolor(edges(e,3)));
            end
        end
    end

    function h=plot_edges_fixed_color_3
        %h=cell(size(edges,1),1);
        if directed
            for e=1:size(edges,1)
                h(e)=arrow3d(xy(edges(e,1:2),1),xy(edges(e,1:2),2),xy(edges(e,1:2),3),...
                    'Color',edgecolor(edges(e,3)),'LineWidth',edgewidth*base_radius*edges(e,3),...
                    'HeadOffset',point_size(edges(e,2))*base_radius*1.01,...
                    'AngularResolution',AngularResolution,...
                    'HeadLength',5/3*options.pointsize(edges(e,2))/edgewidth/base_radius/edges(e,3),...
                    'HeadWidth',options.pointsize(edges(e,2))/edgewidth/base_radius/edges(e,3));
            end
        else
            for e=1:size(edges,1)
                h(e)=plot3(xy(edges(e,1:2),1),xy(edges(e,1:2),2),xy(edges(e,1:2),3),'linewidth',edgewidth,'color',edgecolor(edges(e,3)));
            end
        end
    end



%% set up node colors and node plotting function
if isset(options,'nodecolors')
    map=options.nodecolors;
else
    if size(options.scores,2)==1
        map=([0,0,0;distinguishable_colors(length(unique(options.scores(options.scores~=0))),{'k','w'})]);
    else
        map=[0,0,0;distinguishable_colors(size(options.scores,2),{'k','w'})];
    end
end

if size(map,1)>1
    if size(options.scores,2)==1
        unc=unique(options.scores);
    else
        unc=0:size(options.scores,2);
    end
    switch length(unc)
        case size(map,1)
            nodecolor=@(score) map(unc==score,:);
        case size(map,1)-1
            nodecolor=@(score) map(find(unc==score)+1,:);
        otherwise
            nodecolorlim=options.nodecolorlim;
            if nodecolorlim(1)==nodecolorlim(2)
                nodecolor=@(score) map(1,:);
            else
                nodecolor=@(score) interp1(linspace(nodecolorlim(1),nodecolorlim(2),size(map,1)),map,score,'nearest');
            end
    end
else
    nodecolor=@(score) map;
end

%switch between 2d or 3d plotting
switch size(xy,2)
    case 2
        if size(options.scores,2)==1&&~directed
            plot_node=@(xy,score,shape,pointsize) plot(xy(1),xy(2),shape,'markersize',pointsize,'markerfacecolor',nodecolor(score),'markeredgecolor',nodecolor(score));
        else
            plot_node=@plot_pie;
        end
    case 3
        if size(options.scores,2)>1
            plot_node=@plot_pie_3;
        else
            n_points=AngularResolution;
            [xs,ys,zs]=sphere(n_points);
            plot_node=@(xy,color,shape,pointsize) surf(xy(1)+xs*pointsize,xy(2)+ys*pointsize,xy(3)+zs*pointsize,colorarray(nodecolor(color),n_points+1,n_points+1),'edgecolor','none');
            %plot_node=@(xy,color,shape,pointsize) plot3(xy(1),xy(2),xy(3),shape,'markersize',pointsize,'markerfacecolor',nodecolor(color),'markeredgecolor',nodecolor(color));
        end
    otherwise
        error('need 2 or 3 dimensional coordinates');
end



    function h=plot_pie(xy,score,shape,pointsize)
        
        points=50;
        radius=pointsize*base_radius;
        if sum(score)>0
            shares=score./sum(score);
            last_t=0;
            it_share=find(shares);
            h=zeros(length(it_share),1);
            for it=1:length(it_share)
                end_t=last_t+shares(it_share(it))*points;
                tlist=[last_t, ceil(last_t):floor(end_t), end_t];
                xlist=[0,(radius*cos(tlist*2*pi/points)),0]+xy(1);
                ylist=[0,(radius*sin(tlist*2*pi/points)),0]+xy(2);
                hp=patch(xlist,ylist,nodecolor(it_share(it)),'EdgeColor','none');
                set(hp,'userdata',it_share(it));
                set(get(get(hp,'annotation'),'legendinformation'),'icondisplaystyle','off');
                if options.isset('legendlabels')
                    set(hp,'displayname',options.legendlabels{it_share(it)});
                else
                    set(hp,'displayname',num2str(it_share(it)));
                end
                h(it)=hp;
                last_t=end_t;
            end
        else
            tlist=0:points;
            xlist=xy(1)+radius*cos(tlist*2*pi/points);
            ylist=xy(2)+radius*sin(tlist*2*pi/points);
            hp=patch(xlist,ylist,nodecolor(0),'EdgeColor','none');
            set(hp,'userdata',0)
            set(hp,'displayname','missing')
            set(get(get(hp,'annotation'),'legendinformation'),'icondisplaystyle','off');
            h=hp;
        end
    end

    function h=plot_pie_3(xy,score,shape,pointsize)
        points=50;
        radius=pointsize*base_radius;
        if sum(score)>0
            shares=score./sum(score);
            last_t=0;
            it_share=find(shares);
            h=zeros(length(it_share),1);
            for it=1:length(it_share)
                end_t=last_t+shares(it_share(it))*points;
                tlist=[last_t, ceil(last_t):floor(end_t), end_t];
                xlist=[0,(radius*cos(tlist*2*pi/points)),0]+xy(2);
                ylist=[0,(radius*sin(tlist*2*pi/points)),0]+xy(3);
                hp=patch(xy(1)*ones(size(xlist)),xlist,ylist,nodecolor(it_share(it)),'EdgeColor','none');
                set(hp,'userdata',it_share(it));
                set(get(get(hp,'annotation'),'legendinformation'),'icondisplaystyle','off');
                if options.isset('legendlabels')
                    set(hp,'displayname',options.legendlabels{it_share(it)});
                else
                    set(hp,'displayname',num2str(it_share(it)));
                end
                h(it)=hp;
                last_t=end_t;
            end
        else
            tlist=0:points;
            xlist=xy(2)+radius*cos(tlist*2*pi/points);
            ylist=xy(3)+radius*sin(tlist*2*pi/points);
            hp=patch(xy(1)*ones(size(xlist)),xlist,ylist,nodecolor(0),'EdgeColor','none');
            set(hp,'userdata',0)
            set(hp,'displayname','missing')
            set(get(get(hp,'annotation'),'legendinformation'),'icondisplaystyle','off');
            h=hp;
        end
    end



%% set up axis
if ~is_hold
    cla;
    %     xy_min=min(xy);
    %     xy_max=max(xy);
    %     lims(1:2:2*length(xy_min))=xy_min;
    %     lims(2:2:2*length(xy_max))=xy_max;
    %     axis(1.1*lims);
    caxis(options.edgecolorlim);
    axis off
    axis equal
    %     axis tight
    set(gca,'position',[0.01,0.01,0.98,0.98])
end
hold on

%% plot the edges and nodes
h_edges=plot_edges();

h_nodes=cell(N,1);
scores=options.scores;
for i=1:N
    if point_size(i)
        switch shapes(i,:)
            case ' '
                
                %             case '.'
                %                 h_nodes{i}=plot_node(xy(i,:)-2^-10,scores(i,:),shapes(i,:),point_size(i)*2);
            otherwise
                h_nodes{i}=plot_node(xy(i,:),scores(i,:),shapes(i,:),point_size(i));
        end
    end
end


%% set up legend

% remove edges from legend
edge_annotation=get(h_edges,'Annotation');
if iscell(edge_annotation)
    for i=1:length(edge_annotation)
        set(get(edge_annotation{i},'legendinformation'),'Icondisplaystyle','off');
    end
else
    set(get(edge_annotation,'legendinformation'),'Icondisplaystyle','off');
end

% group nodes with unique score and shape in legend
plotted=find(point_size>0);
if ~isempty(plotted)
    if size(options.scores,2)==1
        [us,~,lind]=unique([options.scores(plotted),double(shapes(plotted))],'rows');
        groups(1:size(us,1))=0;
        for i=1:size(us,1)
            groups(i)=hggroup;
            c_nodes= lind==i;
            set([h_nodes{plotted(c_nodes)}],'parent',groups(i));
        end
        annot=get(groups,'annotation');
        if ~iscell(annot)
            annot={annot};
        end
        for i=1:length(annot)
            if options.isset('legendlabels')
                set(groups(i),'Displayname',options.legendlabels{i});
            else
                set(groups(i),'Displayname',num2str(us(i)));
            end
            set(get(annot{i},'legendinformation'),'icondisplaystyle','on');
        end
    else
        not_in_legend=true(size(options.scores,2)+1,1);
        missing=any(sum(options.scores,2)==0);
        if ~missing
            not_in_legend(1)=false;
        end
        patches=cell(size(options.scores,2)+1,1);
        pos=1;
        while any(not_in_legend)&&pos<=length(plotted)
            i=plotted(pos);
            patch_handles=h_nodes{i};
            if iscell(patch_handles)
                ref_handle=@(j) patch_handles{j};
            else
                ref_handle=@(j) patch_handles(j);
            end
            
            for j=1:length(patch_handles)
                patch_scores=get(ref_handle(j),'userdata');
                if not_in_legend(patch_scores+1)
                    set(get(get(ref_handle(j),'annotation'),'legendinformation'),'icondisplaystyle','on');
                    patches{patch_scores+1}=ref_handle(j);
                    not_in_legend(patch_scores+1)=false;
                end
            end
            pos=pos+1;
        end
        for i=1:length(patches)
            uistack(patches{i},'top')
        end
    end
end

%% clean up and output handles if desired
if ~is_hold
    hold off
end

if nargout
    h_nodes_out=h_nodes;
    h_edges_out=h_edges;
end



end

%% create color array for 3d sphere
function c=colorarray(color,m,n)
c=zeros(m,n,3);
for i=1:m
    for j=1:n
        c(i,j,:)=color;
    end
end
end



