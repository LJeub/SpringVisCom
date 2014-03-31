function rgbmap=colorstr2rgb(colorstring)
% Convert a color specified by colorstring to rgb map
% See LineSpec for colorcodes

% Version:
% Date:
% Author:
% Email:

if ~iscell(colorstring)
    colorstring={colorstring};
end

rgbmap=zeros(length(colorstring),3);
for i=1:length(colorstring)
    rgbmap(i,:)= bitget(find('krgybmcw'==colorstring{i})-1,1:3);
end
