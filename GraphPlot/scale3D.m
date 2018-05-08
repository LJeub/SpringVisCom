function scale_factor=scale3D(axesHandle)

oldAxesUnits=get(axesHandle,'Units');
set(axesHandle,'Units','points');
ax_pos=get(axesHandle,'position');
ax_width=ax_pos(3);
ax_height=ax_pos(4);
set(axesHandle,'Units',oldAxesUnits);
projection=view(axesHandle);
xlim=get(axesHandle,'xlim');
x_length=projection * [xlim(2);0;0;1]- projection*[xlim(1);0;0;1];

vertices = projection * [0, 1, 0, 0, 1, 1, 0, 1;
    0, 0, 1, 0, 1, 0, 1, 1;
    0, 0, 0, 1, 0, 1, 1, 1;
    1, 1, 1, 1, 1, 1, 1, 1];
verticesXY = vertices([1, 2], :);

A = [ 0,  0,  0, -1, +1,  0,  0,  0;
    0,  0, -1,  0,  0, +1,  0,  0;
    0, -1,  0,  0,  0,  0, +1,  0;
    -1,  0,  0,  0,  0,  0,  0, +1];

diagonals = verticesXY * A';

dimensions=max(abs(diagonals),[],2);
if strcmpi(get(axesHandle,'DataAspectRatioMode'),'auto')&&strcmpi(get(axesHandle,'PlotBoxAspectRatioMode'),'auto')
    width=ax_pos(3);
    height=ax_pos(4);
else
    aspectRatio=dimensions(2)*ax_width/(dimensions(1)*ax_height);
    axesAspectRatio=ax_pos(4)/ax_pos(3);
    if aspectRatio>axesAspectRatio
        width = ax_pos(4) / aspectRatio;
        height=ax_pos(4);
    else
        width=ax_pos(3);
        height=ax_pos(3)*aspectRatio;
    end
end

xlen_points = sqrt((x_length(1)/dimensions(1)*width)^2 + (x_length(2)/dimensions(2)*height)^2);
xlen_data = xlim(2)-xlim(1);
scale_factor=xlen_data/xlen_points;
end



