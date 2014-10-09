function a = looseAxis(percent,pa)
% adjust plot axis to be scaled such that there is x% of empty axis around
% the data - only works in 2D
%
% pa can be 'xy', 'x', or 'y' corresponding to which axis will be changed
%
% returns the original axis, not the new axis, do a = axis; to get the new
%   axis

if nargin < 1
    percent = 0.1; pa = 'xy';
elseif nargin < 2
    pa = 'xy';    
end

a = axis;

axis tight;
ax = axis(gca);

ax_d = percent*[ax(2)-ax(1) ax(4)-ax(3)];
if strcmp(pa,'xy')    
    ax_new = [ax(1)-ax_d(1) ax(2)+ax_d(1) ax(3)-ax_d(2) ax(4)+ax_d(2)];
elseif strcmp(pa,'x')
    ax_new = [ax(1)-ax_d(1) ax(2)+ax_d(1) ax(3) ax(4)];
elseif strcmp(pa,'y')
    ax_new = [ax(1) ax(2) ax(3)-ax_d(2) ax(4)+ax_d(2)];
else
    error('Invalid input argument: pa must be ''xy'', ''x'', or ''y''');    
end

axis(ax_new);