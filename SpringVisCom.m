function [xy,xyc,energy,energy_c]=SpringVisCom(A,S,varargin)
% [XY,XYC,ENERGY,ENERGY_C]=SPRINGVISCOM(A,S,OPTIONS)
% calculate coordinates for visualising networks with community structure
%
% Inputs:
%
%   A: adjacency matrix
%
%   S: community assignment:
%           possible forms:
%               group vector gives the community index for each node
%                   [size: (#nodes,1)]
%
%               group matrix: gives the strength of community memberships
%                   for each node to each community, used for overlapping
%                   communities
%                   [size: (#nodes,#groups)]
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
%   xy: array of coordinates for each node [size: (#nodes,dimension)]
%
%   xyc: array of coordinates for each community
%           [size: (#communities,dimension)] (used for e.g. drawForceCPie)
%
%   energy: energy of the spring system in the final layout (same as
%       energy_c if communities_only is set)
%
%   energy_c: energy of the spring system in the community layout (all
%       nodes in a community are located at a single community coordinate
%       given by xyc)
%
% Requires:
%
%   matlab_bgl
%   OptionStruct
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

% Version: 1.2
% Date: Tue 13 May 2014 17:03:19 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk


%% check input
if nargin<2||isempty(S)
    S=ones(size(A,1),1);
end

if ~isequal(size(A,1),size(S,1))
    %transpose group vector if necessary
    if isequal(size(A,1),size(S,2))
        S=S';
    else
        error('group vector has wrong dimension for adjacency matrix')
    end
end


%% option parsing
%set up defaults
options=OptionStruct('epsilon',.01,'c1',0.0001,'c2',0.9,'community_field',1,...
    'dimension',2,'verbose',false,'optimisation_groups',[],...
    'fixed_nodes',[],'progressive',[],'distance_matrix',false,...
    'communities_only',false,'background_field',0);

%parse input
options.set(varargin);

%set local variables for speed
c1=options.c1;
c2=options.c2;
epsilon=options.epsilon;
dim=options.dimension;
fixed_nodes=options.fixed_nodes;
if options.isset('progressive')
    xy=options.progressive;
    if size(fixed_nodes(:),2)==1
        fixed_nodes=[fixed_nodes(:),xy(fixed_nodes,:)];
    end
    %check no two nodes have same coordinates, otherwise add noise
    uxy=unique(xy,'rows');
    while size(uxy,1)<size(xy,1)
        r=.01*randn(size(xy,1),2);
        xy=xy+r;
        uxy=unique(xy,'rows');
    end
else
    xy=rand(length(A),dim);
end

%switch verbose output on or off
if options.verbose
    verbose=@(step,energy,gradient)...
        fprintf('step = %10.5f, energy = %d, grad = %d\n',step,energy,gradient);
else
    verbose=@(step,energy,gradient) [];
end

% parse group input and set up charge matrix
optimisation_groups=optimisation_group_matrix(S);

if size(S,2)==1 %group vector
    %clean group vector
    %nodes coded zero should not belong to any community
    keep=find(S>0);
    [uc,~,S]=unique(S);
    if ismember(0,uc)
        nc=length(uc)-1;
        S=S-1;
    else
        nc=length(uc);
    end
    groups=sparse(keep,S(keep),1,size(S,1),nc);
    charges=charge_matrix(groups);
    
else %group matrix
    [~,groups]=max(S,[],2);
    groups=sparse(1:length(groups),groups,1,size(S,1),size(S,2));
    charges=charge_matrix(S);
end

charges=charges*options.community_field+ones(size(charges))*options.background_field;
for i=1:length(charges)
    charges(i,i)=0;
end

% set up optimisation groups if given
if options.isset('optimisation_groups')
    optimisation_groups=optimisation_group_matrix(options.optimisation_groups);
end

% set up fixed nodes
if options.isset('fixed_nodes')
    optimisation_groups(fixed_nodes(:,1),:)=0;
end
optimisation_groups_sizes=sum(optimisation_groups,1);
opt_nc=size(optimisation_groups,2);

% set up spring system
if options.distance_matrix
    distances=A;
else
    distances=distance_matrix(A);
end
sp_sys=spring_system(distances,charges,groups,fixed_nodes);


%% solve the community layout problem

% initialise community coordinates
if ~options.isset('progressive')
    xyc=rand(sp_sys.nc,dim);
elseif options.communities_only
    xyc=options.progressive;
else
    xyc=[];
end

% optimise community layout (unless progressive is specified and layout is
% not communities_only)
if ~isempty(xyc)
    L=eye(dim*sp_sys.nc);
    step=1;
    gg=@(xy) sp_sys.gradient_com(xy);
    ee=@(xy) sp_sys.energy_com(xy);
    grad=gg(xyc);
    while max(norm2(grad))>epsilon
        e_old=ee(xyc);
        try
            [xyc,grad,step,L]=BFGS_update(xyc,grad,1,L,c1,c2,gg,ee);
            verbose(step,e_old-ee(xyc),max(abs(grad(:))));
        catch err
            switch err.message
                case 'stalled'
                    if options.isset('verbose')
                        warning('no further improvement found')
                    end
                    break
                otherwise
                    rethrow(err);
            end
        end
    end
    
    energy_c=ee(xyc);
    if options.communities_only
        xy=zeros(sp_sys.n,dim);
        for i=1:sp_sys.nc
            ind=find(sp_sys.communities(:,i));
            xy(ind,:)=repmat(xyc(i,:),length(ind),1);
        end
        energy=energy_c;
        return;
    else
        % convert to coordinates for all nodes with perturbation
        xy=zeros(sp_sys.n,dim);
        for i=1:sp_sys.nc
            ind=find(sp_sys.communities(:,i));
            xy(ind,:)=repmat(xyc(i,:),length(ind),1)+randn(length(ind),dim)*.01;
        end
    end
end


%% set coordinates of fixed nodes
if options.isset('fixed_nodes')
    xy(fixed_nodes(:,1),:)=fixed_nodes(:,2:end);
end


%% solve node layout problem
clear L;
sp_sys.set_xy(xy);
delta=norm2(sp_sys.gradient);
delta_com=(delta'*optimisation_groups)./optimisation_groups_sizes;

for i=1:opt_nc
    L{i}=eye(dim*optimisation_groups_sizes(i));
end
[~,com]=max(delta_com);
%com=randsample(sp_sys.nc,1,true,delta_com);
nodes=find(optimisation_groups(:,com));
step=ones(length(L),1);

it=1;
while max(delta)>epsilon
    grad=sp_sys.gradient([],nodes);
    e_old=sp_sys.energy;
    
    it=it+1;
    xy_nodes=sp_sys.xy(nodes,:);
    gg=@(xy) sp_sys.gradient(xy,nodes);
    ee=@(xy) sp_sys.node_energy(xy,nodes);
    try
        [xy_nodes,grad,step(com),L{com}]=BFGS_update(xy_nodes,grad,1,L{com},c1,c2,gg,ee);
        sp_sys.move(nodes,xy_nodes,grad);
    catch err
        switch err.message
            case 'stalled'
                warning('no further improvement found')
                break
            otherwise
                rethrow(err);
        end
    end
    
    delta=norm2(sp_sys.gradient);
    delta_com=(delta'*optimisation_groups)./optimisation_groups_sizes;
    [~,com]=max(delta_com);
    %com=randsample(sp_sys.nc,1,true,delta_com);
    nodes=find(optimisation_groups(:,com));
    verbose(step(com),e_old-sp_sys.energy,max(delta));
end

%% output node coordinates
xy=sp_sys.xy;
energy=sp_sys.energy;
end

%% element-wise 2-norm
function norm=norm2(v)
norm=sqrt(sum(v.^2,2));
end

%% compute charge matrix
function charges=charge_matrix(groups)
charges=1-div_0(groups*groups',norm2(groups)*norm2(groups)');
end

%% create optimisation group matrix
function groups=optimisation_group_matrix(groups)
if size(groups,1)==1||size(groups,2)==1
    [~,~,groups]=unique(groups);
    groups=sparse(1:length(groups),groups,1);
else
    groups=double(logical(groups));
    % ensure all nodes appear in at least one optimisation group
    missing=find(sum(groups,2)==0);
    groups=[groups,sparse(missing,1:length(missing),1)];
end
end

%% convert adjacency matrix to distance matrix
function distances=distance_matrix(A)
[ii,jj,vv]=find(A);
distances=all_shortest_paths(sparse(ii,jj,(max(vv)-vv+1),size(A,1),size(A,2)));
end


