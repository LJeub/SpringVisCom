function MultilayerGraphPlot(xy,A,varargin)



layer_width=length(A{1});
n_layers=length(A);
N=layer_width*n_layers;
for i=1:n_layers
    xyl((1:layer_width)+(i-1)*layer_width,1)=i;
    xyl((1:layer_width)+(i-1)*layer_width,2:3)=xy;
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

options=OptionStruct('alpha',1,'scores',ones(N,1),'shapes','.','pointsize',7,'edgewidth',1,'nodecolors',[],'nodecolorlim',[],...
    'edgecolors',[],'edgecolorlim',[],'randedges',0,'edgethreshold',[],'legendlabels',[],'axispadding',0.1,'layercolor',[0.5,0.5,0.5],'layerlabels',[],'layeralpha',0.1,'aspectratio',[1,1,1],'view',[-15,30]);
options=options.set(varargin);
options.pointsize=options.pointsize*double(logical(k));

if iscell(options.scores)
    scores=[];
    for i=1:length(options.scores)
        scores=[scores;options.scores{i}];
    end
    options.scores=scores;
else
    options.scores=options.scores(:);
    options.scores=sparse(1:length(options.scores),options.scores,1);
end

plot_options=options.struct;
plot_options=rmfield(plot_options,{'axispadding','layercolor','layeralpha','layerlabels','aspectratio','view'});
GraphPlot(xyl,AS,plot_options);



view(options.view);
xmin=1-options.axispadding;
xmax=n_layers+options.axispadding;
ymin=min(xy(:,1));
ymin=ymin-options.axispadding;
ymax=max(xy(:,1));
ymax=ymax+options.axispadding;
zmin=min(xy(:,2));
zmin=zmin-options.axispadding;
zmax=max(xy(:,2));
zmax=zmax+options.axispadding;

axis([xmin,xmax,ymin,ymax,zmin,zmax]);

for i=1:n_layers
    patch_h(i)=patch(2^-10+[i,i,i,i],[ymin,ymax,ymax,ymin],[zmin,zmin,zmax,zmax],options.layercolor);
end
set(patch_h,'facealpha',options.layeralpha);
set(gca,'dataaspectratio',options.aspectratio)

if options.isset('layerlabels')
    set(gca,'xtick',xmin:xmax)
    set(gca,'xticklabel',options.layerlabels);
    set(gca,'xticklabelrotation',-45)
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    axis on
end

    