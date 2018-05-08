function MultilayerGraphPlot(xy,A,varargin)
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
GraphPlot(xyl,AS,plot_options);

if options.drawlayers
    for i=1:n_layers
        patch_h(i)=patch([i,i,i,i]/aspect_ratio(1),[ymin,ymax,ymax,ymin],...
            [zmin,zmin,zmax,zmax],options.layercolor,'linestyle','none');
    end
    set(patch_h,'facealpha',options.layeralpha);
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


