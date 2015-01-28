function xy=MultilayerSpringVisCom(A,S,varargin)



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
    
    D(ind,ind)=distance_matrix(A{i});
    C(ind,ind)=charge_matrix(groups{i});
    
    node_ind(ind)=1:layer_width;
end

D(isnan(D))=0;
D(D==inf)=0;

[~,xy]=SpringVisCom(D,node_ind,'distance_matrix',true,'verbose',true,'communities_only',true,'charge_matrix',C);


    function distances=distance_matrix(A)
        [ii,jj,vv]=find(A);
        distances=all_shortest_paths(sparse(ii,jj,(max(vv)-vv+1),size(A,1),size(A,2)));
    end

    function charges=charge_matrix(groups)
        charges=1-(groups*groups')./(norm2(groups)*norm2(groups)'+eps);
    end

    function norm=norm2(v)
        norm=sqrt(sum(v.^2,2));
    end
end