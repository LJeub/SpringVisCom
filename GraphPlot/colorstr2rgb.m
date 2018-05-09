function rgbmap=colorstr2rgb(colorstring)
% Convert a color specified by colorstring to rgb map
% See LineSpec for colorcodes

% Version: 1.3
% Date: Wed  9 May 2018 15:25:32 CEST
% Author: Lucas Jeub
% Email: ljeub@iu.edu

if ~iscell(colorstring)
    colorstring={colorstring};
end

rgbmap=zeros(length(colorstring),3);
for i=1:length(colorstring)
    rgbmap(i,:)= bitget(find('krgybmcw'==colorstring{i})-1,1:3);
end
