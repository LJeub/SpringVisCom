classdef spring_system < matlab.mixin.Copyable
% SPRING_SYSTEM spring system class for SpringCoordC

% Version: 1.2
% Date: Tue 13 May 2014 17:03:19 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk
    
    %% properties
    properties (SetAccess=private)
        n;
        nc;
        nf;
        dim;
        
        communities;
        sizes_com;
        fixed_nodes;
        
        xy;
        xyc;
        xyf;
    end
    
    %% private properties
    properties (Access=private)
        
        k_com;
        kl_com;
        kl2_com;
        charges_com;
       
        charges;
        spring_length;
        spring_strength;
    end
    
    properties (Transient=true,Access=private)
        distances;
        distances_com;
        curr_gradient;
        curr_energy;
        curr_gradient_com;
        curr_energy_com;
        
        K;
        K_com;
    end
    
    %% methods
    methods
        
        %% constructor
        function obj=spring_system(distances,charges,groups,fixed_nodes)
            if nargin==4
                % distances: distances between nodes in network
                % charges: matrix of charges between pairs of nodes
                % groups: groups of nodes for optimisation
                % fixed_nodes: nodes to keep fixed at given coordinates
                obj.n=length(distances);
                
                obj.spring_length=distances;
                obj.spring_strength=1./distances.^2;
                obj.spring_strength(distances==0)=0;
                obj.spring_length(obj.spring_length==inf)=0;
                
                
                % nodes that are fixed
                if ~isempty(fixed_nodes)
                    obj.fixed_nodes=fixed_nodes(:,1);
                    obj.nf=size(fixed_nodes,1);
                    obj.xyf=fixed_nodes(:,2:end);
                else
                    obj.fixed_nodes=[];
                    obj.nf=0;
                    obj.xyf=[];
                end
                
                
                % number of groups
                obj.nc=size(groups,2);
                
                % fixed nodes get split out into their own groups
                if ~isempty(fixed_nodes)
                    groups(fixed_nodes(:,1),:)=0;
                    g_f=sparse(fixed_nodes(:,1),1:obj.nf,true,obj.n,obj.nf);
                    groups=[groups,g_f];
                end
                obj.communities=groups;
                obj.sizes_com=sum(groups,1);
                
                %set up charges
                obj.charges=charges;
                
                %set up coefficients for calculating community gradient and energy
                obj.k_com=obj.communities'*obj.spring_strength...
                    *obj.communities;
                obj.kl_com=obj.communities'*(obj.spring_strength...
                    .*obj.spring_length)*obj.communities;
                obj.kl2_com=obj.spring_strength.*(obj.spring_length.^2);
                for i=1:obj.n
                    obj.kl2_com(i,i)=0;
                end
                obj.kl2_com=sum(obj.kl2_com(:));
                
                obj.charges_com=obj.communities'*obj.charges*obj.communities;
                
            else
                error('need 4 input arguments to set up spring system')
            end
        end
        
        %% set coordinates
        function set_xy(obj,xy)
            if size(xy,1)~=obj.n
                error('coordinates are the wrong size');
            end
            obj.xy=xy;
            if isempty(obj.dim)
                obj.dim=size(xy,2);
            elseif ~obj.dim==size(xy,2)
                error('coordinates have the wrong dimension');
            end
            obj.distances=squareform(pdist(xy));
            obj.K=obj.spring_strength-(obj.spring_length...
                .*obj.spring_strength)./obj.distances...
                -obj.charges./obj.distances.^3;
            for i=1:obj.n
                obj.K(i,i)=0;
            end
            obj.curr_gradient=obj.gradient(xy);
            obj.curr_energy=obj.energy(xy);
        end
        
        
        function set_xyc(obj,xyc)
            if size(xyc,1)~=obj.nc
                error('coordinates are the wrong size');
            end
            obj.xyc=xyc;
            if isempty(obj.dim)
                obj.dim=size(xyc,2);
            elseif ~obj.dim==size(xyc,2)
                error('coordinates have the wrong dimension')
            end
            obj.distances_com=squareform(pdist(xyc));
            obj.K_com=obj.k_com-obj.kl_com./obj.distances_com...
                -obj.charges_com./obj.distances_com.^3;
            for i=1:obj.nc
                obj.K_com(i,i)=0;
            end
            obj.curr_gradient_com=obj.gradient_com(xyc);
            obj.curr_energy_com=obj.energy_com(xyc);
        end
        
        
        
        %% compute gradient
        function grad=gradient(obj,xy,i)
            if nargin <3
                i=1:obj.n;
            end
            
            if nargin < 2 ||isempty(xy)
                if isempty(obj.xy)
                    error('spring system does not have associated coordinates')
                else
                    grad=obj.curr_gradient(i,:);
                end
            else
                if size(xy,1)<obj.n
                    % new coordinates given only for moved nodes
                    if isempty(obj.xy)
                        error('spring system does not have associated coordinates')
                    end
                    xy_all=obj.xy;
                    xy_all(i,:)=xy;
                    xy=xy_all;
                end
                grad=zeros(length(i),size(xy,2));
                %dist=sqrt(sqdistance(xy',xy(i,:)'));
                dist=pdist2(xy,xy(i,:));
                
                K=obj.spring_strength(:,i)-(obj.spring_length(:,i)...
                    .*obj.spring_strength(:,i))./dist...
                    -obj.charges(:,i)./dist.^3;
                for j=1:length(i)
                    K(i(j),j)=0;
                end
                
                for j=1:length(i)
                    for k=1:size(xy,2)
                        grad(j,k)=sum((xy(i(j),k)-xy(:,k)).*K(:,j));
                    end
                end
            end
        end
        
        function grad=gradient_com(obj,xyc,i)
            if nargin<3
                i=1:obj.nc;
            end
            if nargin<2||isempty(xyc)
                if isempty(obj.xyc)
                    error('spring system does not have associated community coordinates');
                else
                    grad=obj.curr_gradient_com(i,:);
                end
            else
                grad=zeros(length(i),size(xyc,2));
                xyc_all=[xyc;obj.xyf];
                %dist=sqrt(sqdistance(xyc',xyc(i,:)'));
                dist=pdist2(xyc_all,xyc(i,:));
                K_com=obj.k_com(:,i)-obj.kl_com(:,i)./dist-obj.charges_com(:,i)./dist.^3;
                for j=1:length(i)
                    K_com(i(j),j)=0;
                end
                for j=1:length(i)
                    for k=1:size(xyc,2)
                        grad(j,k)=sum((xyc(i(j),k)-xyc_all(:,k)).*K_com(:,j));
                    end
                end
            end
        end
        
        %% compute energy
        function e=energy(obj,xy,nodes)
            if nargin<2
                if isempty(obj.xy)
                    error('spring system does not have associated coordinates')
                else
                    e=obj.curr_energy;
                end
            else
                if nargin>2
                    if isempty(obj.xy)
                        error('spring system does not have associated coordinates')
                    else
                        xy_all=obj.xy;
                        xy_all(nodes,:)=xy;
                        d=pdist2(xy_all,xy);
                        dist=obj.distances;
                        dist(nodes,:)=d';
                        dist(:,nodes)=d;
                    end
                else
                    dist=squareform(pdist(xy));
                end
                e=0.25*obj.spring_strength.*(dist-obj.spring_length).^2 ...
                    +0.5*obj.charges./dist;
                for i=1:obj.n
                    e(i,i)=0;
                end
                e=sum(e(:));
            end
        end
        
        function e=energy_com(obj,xy)
            if nargin<2
                if isempty(obj.xyc)
                    error('no community coordinates')
                else
                    e=obj.curr_energy_com;
                end
            else
                xy=[xy;obj.xyf];
                dist=squareform(pdist(xy));
                e=0.25*dist.^2.*obj.k_com-0.5.*dist.*obj.kl_com...
                    +0.5*obj.charges_com./dist;
                for i=1:obj.nc+obj.nf
                    e(i,i)=0;
                end
                e=sum(e(:))+0.25*obj.kl2_com;
            end
        end
        
        function e=node_energy(obj,xy,nodes)
            xy_all=obj.xy;
            xy_all(nodes,:)=xy;
            d=pdist2(xy_all,xy);
            e=0.5*obj.spring_strength(:,nodes).*(d-obj.spring_length(:,nodes)).^2 ...
                +(obj.charges(:,nodes)./d);
            for i=1:length(nodes)
                e(nodes(1:i),i)=0;
            end
            e=sum(e(:));
        end
        
        
        %% move node and update gradient, energy
        function move(obj,nodes, xy, grad)
            % move node node to new coordinates xy
            if size(xy,1)~=length(nodes)
                error('coordinates have wrong size');
            end
            
            if size(xy,2)~=obj.dim
                error('coordinates have wrong dimension');
            end
            
            if nargin<3
                grad=obj.gradient(xy,nodes);
            end
            
            xyold=obj.xy(nodes,:);
            
            % substract old contributions of moved nodes to gradients
            for i=1:length(nodes)
                for j=1:obj.dim
                    obj.curr_gradient(:,j)=obj.curr_gradient(:,j)-(obj.xy(:,j)-xyold(i,j)).*obj.K(:,nodes(i));
                end
            end
            
            % substract old contributions of moved nodes to energy
            delta_e=0.5*obj.spring_strength(:,nodes).*(obj.distances(:,nodes)-obj.spring_length(:,nodes)).^2 ...
                +(obj.charges(:,nodes)./obj.distances(:,nodes));
            
            for i=1:length(nodes)
                delta_e(nodes(1:i),i)=0;
            end
            
            obj.curr_energy=obj.curr_energy-sum(delta_e(:));
            
            
            % update coordinates
            obj.xy(nodes,:)=xy;
            % changed distances
            d=pdist2(obj.xy,xy);% update distances
            obj.distances(nodes,:)=d';
            obj.distances(:,nodes)=d;
            
            % update contributions
            k=obj.spring_strength(:,nodes)-(obj.spring_length(:,nodes)...
                .*obj.spring_strength(:,nodes))./d...
                -(obj.charges(:,nodes)./d.^3);
            for i=1:length(nodes)
                k(nodes(i),i)=0;
            end
            obj.K(:,nodes)=k;
            obj.K(nodes,:)=k';
            
            % add new contributions to gradient
            for i=1:length(nodes)
                for j=1:obj.dim
                    obj.curr_gradient(:,j)=obj.curr_gradient(:,j)+(obj.xy(:,j)-xy(i,j)).*k(:,i);
                end
            end
            
            
            % calculate gradient of moved nodes
            if nargin<4
                obj.curr_gradient(nodes,:)=obj.gradient(xy,nodes);
            else
                obj.curr_gradient(nodes,:)=grad;
            end
            
            % update energy
            delta_e=0.5*obj.spring_strength(:,nodes).*(d-obj.spring_length(:,nodes)).^2 ...
                +(obj.charges(:,nodes)./d);
            for i=1:length(nodes)
                delta_e(nodes(1:i),i)=0;
            end
            obj.curr_energy=obj.curr_energy+sum(delta_e(:));
        end
    end
end
