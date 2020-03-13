function contourmode(x,y,mode,dB,xyrange);

% Produces a contour plot (in dB) of one field component of the
% mode of an optical waveguide.
% 
% USAGE:
% 
% contourmode(x,y,mode);
% contourmode(x,y,mode,dBrange);
% contourmode(x,y,mode,dBrange,xyrange);
% 
% INPUT:
% 
% x,y - vectors describing horizontal and vertical grid points
% mode - the mode or field component to be plotted
% dBrange - contour levels to plot (in dB), with 0 dB corresponding
%   to the level |mode| = 1. default = (0:-3:-45)
% xyrange - axis range to use (optional)
%
% EXAMPLE:  Make a contour plot of the magnetic field component Hx,
% with contours from 0 dB down to -50 dB, relative to the maximum
% value, in 5 dB increments. 
%
%     contourmode(x,y,Hx/max(abs(Hx(:))),(0:-5:-50));
%
% NOTES:  
%
% (1) This function uses the current color map to determine the
%     colors of each contour, with 0 dB corresponding to the
%     maximum color and -dbmax corresponding to the minimum color.
%     You can use the 'colormap' command to change the current
%     color map.
% (2) The aspect ratio of the plot box is automatically adjusted so
%     that the horizontal and vertical scales are equal.
% (3) The mode is not normalized or scaled in any way.

x = real(x);
y = real(y);

if (nargin < 5)
  xyrange = [min(x),max(x),min(y),max(y)];
end

if (size(mode) == [length(x)-1,length(y)-1])
    x = (x(1:end-1) + x(2:end))/2;
    y = (y(1:end-1) + y(2:end))/2;
end

if (nargin < 4) || isempty(dB)
  dB = (0:-3:-45);
end

% Compute and plot contours
surf(x,y,20*log10(abs(transpose(mode))));

axis(xyrange);
set(gca,'PlotBoxAspectRatio',[xyrange(2)-xyrange(1) xyrange(4)-xyrange(3) 1],...
        'Box','on');
