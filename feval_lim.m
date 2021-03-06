function YY = feval_lim(xx, yy, xx_lim, yy_lim)
%FEVAL_LIM Restrict values to be monotonic relative to data
%
% YY = FEVAL_LIM(xx,yy,xx_lim,yy_lim) takes data that has been evaluated at
% (xx,yy), and ensures that it is monotonic with respect to 
% (xx_lim, yy_lim).
%
% Examples:
% This illustrates the limiting of data (XX,YY) to the data (x_lim, y_lim).
% The limited output is plotted in blue and the (x_lim, y_lim) are plotted
% as red crosses. Note that the data in [1,2] and [4,5] are now strictly
% monotonic.
%     x_lim = 0:6;
%     y_lim = sin(x_lim);
%     XX = linspace(0,6);
%     YY = sin(XX);
%     Y_monotonic = feval_lim(XX,YY,x_lim,y_lim);
%     plot(XX,Y_monotonic,'b-',x_lim,y_lim,'r+')

% Now evaluate the y's where they're supposed to be!
xx_N = length(xx);

YY = yy(:);

BREAKS = xx_lim;
COEFS = yy_lim(:);

[~,BIN] = histc(xx,BREAKS);
yy_max = [max(COEFS(1:end-1),COEFS(2:end));COEFS(end)];
yy_min = [min(COEFS(1:end-1),COEFS(2:end));COEFS(end)];
YY = max(min(YY,yy_max(BIN)),yy_min(BIN));

YY = reshape(YY,size(yy));


end % function
