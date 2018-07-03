function xy=MultilayerSpringVisCom(A,S,varargin)
% [XY]=MULTILAYERSPRINGVISCOM(A,S,OPTIONS)
%
% calculate coordinates for visualising multilayer networks with community
% structure
%
% Inputs:
%
%   A: cell array of adjacency matrices for each layer
%
%   S: community assignment:
%           possible forms:
%               group vector matrix: each column gives the community index 
%                   for each node in a layer
%                   [size: (#nodes,#layers)]
%
%               group matrix array: cell array of group matrices for each 
%                   layer. Each matrix gives the strength of community 
%                   memberships for each node to each community, used for 
%                   overlapping communities
%                   [size: {(#nodes,#groups)}x#layers]
%
%               emtpty: ([]) results in the usual Kamada-Kawai layout
%                   without community structure
%
%
%   OPTIONS:
%           epsilon: [0.01] tolerance for termination criterion
%
%           community_field: [1] repulsive charge between nodes in different
%               communities, increase to separate communities more
%
%           background_field: [0] repulsive charge between all nodes (useful
%               in combination with distance_matrix=true, when the distance
%               matrix is sparse)
%
%           dimension: [2] dimension of the layout space
%
%           distance_matrix: [false] set to true to treat A as a distance
%               matrix rather than an adjacency matrix
%
%           verbose: [false] set to true to display improvement at each step
%
%           optimisation_groups: [S] groups used for optimisation
%               (group vector or group matrix)
%
%           fixed_nodes: [] format [node, x, y,...] keep these nodes fixed
%               at given coordinates
%
%           progressive: [] initial coordinates for all nodes, when
%               progressive is specified, fixed_nodes can be given as list
%               of node indices, coordinates specified in fixed_nodes take
%               precendence over those in progressive, if communities_only
%               is set, give coordinates for each community
%
%           communities_only: [false] set to true to compute only community
%               coordinates
%
%           c1: [0.0001] first constant in Wolfe-conditions for line-search
%
%           c2: [0.9] second constant in Wolfe-conditions for line-search
%
%   options are specified as a comma-separated list of key-value pairs or
%   as a struct of options
%
% Outputs:
%
%   xy: array of coordinates for each physical node [size: (#nodes,dimension)]
%
% Requires:
%
%   matlab_bgl
%   spring_system
%   BFGS_update
%   wolfe_step
%   div_0
%
% The node layout is computed by minimising an energy function which is
% based on a spring system (see Tomihisa Kamada and Satoru Kawai:
% "An Algorithm For Drawing General Undirected Graphs"), where nodes are
% charged particles, such that nodes in different communities repel each
% other.

% Version: 1.3.1
% Date: Tue  3 Jul 2018 12:38:16 CEST
% Author: Lucas Jeub
% Email: lucasjeub@gmail.com
parseArgs=inputParser();
parseArgs.KeepUnmatched=true;
addParameter(parseArgs,'background_field',0);
addParameter(parseArgs,'community_field',1);
addParameter(parseArgs,'distance_matrix',false);
parse(parseArgs,varargin{:});

options=parseArgs.Results;
other_options=parseArgs.Unmatched;
n_layers=length(A);
layer_width=length(A{1});

if nargin<2||isempty(S)
    S=cell(n_layers);
    for i=1:n_layers
        S{i}=sparse(layer_width,1);
    end
end

if ~iscell(S)
    for i=1:n_layers
        [~,~,e]=unique(S(:,i));
        groups{i}=sparse(1:layer_width,e,1,layer_width,max(e));
    end
else
    groups=S;
end

for i=1:n_layers
    ind=(1:layer_width)+(i-1)*layer_width;
    if options.distance_matrix
        D(ind,ind)=A{i};
    else
        D(ind,ind)=distance_matrix(A{i});
    end
    C(ind,ind)=charge_matrix(groups{i},A{i});
    node_ind(ind)=1:layer_width;
end

D(isnan(D))=100;
D(D==inf)=100;
D=(D+D')/2;

[~,xy]=SpringVisCom(D,node_ind,'distance_matrix',true,...
    'communities_only',true,'charge_matrix',C, other_options);


    function distances=distance_matrix(A)
        [ii,jj,vv]=find(A);
        distances=all_shortest_paths(sparse(ii,jj,(max(vv)-vv+1),size(A,1),size(A,2)));
    end

    function charges=charge_matrix(groups, A)
        charges=1-(groups*groups')./(norm2(groups)*norm2(groups)'+eps);
        background=ones(size(A));
        background(A>0)=0;
        charges=charges.*(options.community_field(:)*options.community_field(:)')+...
        options.background_field(:)*options.background_field(:)'.*background;
    end

    function norm=norm2(v)
        norm=sqrt(sum(v.^2,2));
    end
end
