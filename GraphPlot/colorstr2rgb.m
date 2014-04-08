function rgbmap=colorstr2rgb(colorstring)
% Convert a color specified by colorstring to rgb map
% See LineSpec for colorcodes

% Version: 1.1
% Date: Mon 31 Mar 2014 17:39:40 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

if ~iscell(colorstring)
    colorstring={colorstring};
end

rgbmap=zeros(length(colorstring),3);
for i=1:length(colorstring)
    rgbmap(i,:)= bitget(find('krgybmcw'==colorstring{i})-1,1:3);
end
