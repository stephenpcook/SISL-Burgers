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

YY = yy;

BREAKS = xx_lim;
COEFS = yy_lim;

for jj = 1:xx_N
    % Find the correct interval. Could also be done with histc?
    interval = 1;
    while xx(jj)>BREAKS(interval + 1)
        interval = interval + 1;
    end % while
    y_l = COEFS(interval);
    y_r = COEFS(interval+1);
    y_max = max(y_l,y_r);
    y_min = min(y_l,y_r);
    if YY(jj) > y_max
        YY(jj) = y_max;
        % Warning messages if appropriate verbatim param
        %warning(['Overshoot limited in interval (',num2str(x(interval)),...
        %    ', ',num2str(x(interval+1)),')'])
    elseif YY(jj) < y_min
        YY(jj) = y_min;
        %warning(['Undershoot limited in interval (',num2str(x(interval)),...
        %   ', ',num2str(x(interval+1)),')'])
        % Warning messages if appropriate verbatim param
    end % if
end % for jj=1:y_N
end % function
