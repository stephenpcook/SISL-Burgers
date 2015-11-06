function YY = feval_lim(xx, yy, xx_lim, yy_lim)
% Takes a function that has been evaluated at (xx,yy), and ensures that it
% is monotonic with respect to (xx_lim, yy_lim)
%
% YY = feval_lim(xx, yy, xx_lim, yy_lim)

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
