function MultilayerGraphPlot(xy,A,varargin)
layer_width=length(A{1});
n_layers=length(A);
N=layer_width*n_layers;

options=OptionStruct('alpha',1,'scores',ones(N,1),'shapes','.','pointsize',7,'edgewidth',1,'nodecolors',[],'nodecolorlim',[],...
    'edgecolors',[],'edgecolorlim',[],'randedges',0,'edgethreshold',[],'legendlabels',[],'axispadding',0.1,'layercolor',[0.5,0.5,0.5],'layerlabels',[],'labelrotation',30,'labelfont','Times','layeralpha',0.1,'aspectratio',[1,1,1],'view',[-15,30],'baseradius',1);
options=options.set(varargin);

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


options.pointsize=options.pointsize(:);
%options.pointsize=options.pointsize*double(logical(k));

if iscell(options.scores)
    scores=[];
    for i=1:length(options.scores)
        scores=[scores;options.scores{i}];
    end
    options.scores=scores;
else
    options.scores=options.scores(:);
    if length(options.scores)<N
        options.scores=sparse(options.scores,1,1,N,1);
    end
end

plot_options=options.struct;
plot_options=rmfield(plot_options,{'axispadding','layercolor','layeralpha','layerlabels','labelrotation','aspectratio','view','labelfont'});
GraphPlot(xyl,AS,plot_options);



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



for i=1:n_layers
    patch_h(i)=patch([i,i,i,i]/aspect_ratio(1),[ymin,ymax,ymax,ymin],[zmin,zmin,zmax,zmax],options.layercolor,'linestyle','none');
end
set(patch_h,'facealpha',options.layeralpha);
axis equal
axis([xmin,xmax,ymin,ymax,zmin,zmax]);
set(gca,'clipping','off')
if options.isset('layerlabels')
%     set(gca,'xtick',xmin:xmax)
%     set(gca,'xticklabel',options.layerlabels);
%     set(gca,'xticklabelrotation',-45)
%     set(gca,'ytick',[])
%     set(gca,'ztick',[])
%     axis on

    for i=1:length(options.layerlabels)
        text(i/aspect_ratio(1),ymax,zmin,options.layerlabels{i},....
            'HorizontalAlignment','left','VerticalAlignment','cap',...
            'rotation',options.labelrotation,'Fontname',options.labelfont);
    end
end

set(gca,'ydir','reverse')

    
