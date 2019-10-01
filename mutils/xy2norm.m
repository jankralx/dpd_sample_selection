%author: Girish Ratanpal, ECE, UVa. 
%transform axes units to normalized units for 2-D figures only. 
%works for linear,log scales and also reverse axes.
%DATE: JAN. 6, 2006. previous version: aug 12, 2005.

% Kral 2019-06-03: function renamed from y_to_norm_v2 to xy2norm
% first parameter is 'x' or 'y' which determines the axis

%FUNCTION DESCRIPTION:
% function [y_norm] = xy2norm('y',y_begin,y_end)
% function returns a 1x2 array, y_norm, with normalized units y_begin and
% y_end for the annotation object. 
%     
% function arguments:
%1. axis: 'x' or 'y' selects axis
%2. xy_begin: enter the point where the annotation object begins on the axis 
%         using axis units
%3. xy_end: end of annotation object, again in axis units. 
%         

%EXAMPLE: first plot the graph on the figure. 
%then use the following command for placing an arrow:
%h_object =
%annotation('arrow',xy2norm('x',x_begin,x_end),xy2norm('y',y_begin,y_end));

function [xy_norm] = xy2norm(axis, xy_begin, xy_end)

if nargin ~= 3
     error('Wrong number of input arguments! xy2norm(axis,xy_begin,xy_end)')
end

if axis ~= 'x'
    axis = 'y';
end

h_axes = get(gcf,'CurrentAxes');    %get axes handle.
axesoffsets = get(h_axes,'Position');%get axes position on the figure. 
xy_axislimits = get(gca, sprintf('%slim',axis));     %get axes extremeties.
xy_dir = get(gca,sprintf('%sdir',axis));
xy_scale = get(gca,sprintf('%sScale',upper(axis)));

%get axes length
xy_axislength_lin = abs(xy_axislimits(2) - xy_axislimits(1));


if strcmp(xy_dir,'normal')      %axis not reversed
  if strcmp(xy_scale,'log')
    %get axes length in log scale.
    xy_axislength_log = abs(log10(xy_axislimits(2)) - log10(xy_axislimits(1)));  
    
    %normalized distance from the lower left corner of figure. 
    xy_begin_norm = axesoffsets(2)+axesoffsets(4)*abs(log10(xy_begin)-log10(xy_axislimits(1)))/(xy_axislength_log);
    xy_end_norm = axesoffsets(2)+axesoffsets(4)*abs(log10(xy_end)-log10(xy_axislimits(1)))/(xy_axislength_log);

    xy_norm = [xy_begin_norm xy_end_norm];
  elseif strcmp(xy_scale,'linear')%linear scale.
    %normalized distance from the lower left corner of figure.  
    xy_begin_norm = axesoffsets(2)+axesoffsets(4)*abs((xy_begin-xy_axislimits(1))/xy_axislength_lin);
    xy_end_norm = axesoffsets(2)+axesoffsets(4)*abs((xy_end-xy_axislimits(1))/xy_axislength_lin);

    xy_norm = [xy_begin_norm xy_end_norm];  
  else
      error('use only lin or log in quotes for scale')
  end   
   
elseif strcmp(ydir,'reverse')  %axis is reversed
    if strcmp(xy_scale,'log')
        %get axes length in log scale.
        xy_axislength_log = abs(log10(xy_axislimits(2)) - log10(xy_axislimits(1)));
        %normalized distance from the lower left corner of figure. 
        xy_begin_norm = axesoffsets(2)+axesoffsets(4)*abs(log10(xy_axislimits(2))-log10(xy_begin))/(xy_axislength_log);
        xy_end_norm = axesoffsets(2)+axesoffsets(4)*abs(log10(xy_axislimits(2))-log10(xy_end))/(xy_axislength_log);

        xy_norm = [xy_begin_norm xy_end_norm]; 
    elseif strcmp(xy_scale,'linear')
        %normalized distance from the lower left corner of figure. 
        xy_begin_norm = axesoffsets(2)+axesoffsets(4)*abs((xy_axislimits(2)-xy_begin)/xy_axislength_lin);
        xy_end_norm = axesoffsets(2)+axesoffsets(4)*abs((xy_axislimits(2)-xy_end)/xy_axislength_lin);

        xy_norm = [xy_begin_norm xy_end_norm];
    else
        error('use only lin or log in quotes for scale')
    end
else
    error('use only r or nr in quotes for reverse')
end