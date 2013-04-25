% Creates a plot for use in reports

clf
clear all
x = [0,0.3,0.5,0.7,1];
M = [0.5,0.3,0.25,0.25,0.4];
M=M./trapz(x,M);
plot(x,M)
ylim([0,2])
Monitor.type='points';
Monitor.x = x;
Monitor.M = M;
hold on
xi=Eqd1dExact(x,Monitor);
ML=interp1(x,M,xi);
for ii = 1:5
  plot([xi(ii) xi(ii)],[ML(ii) 0],'r-')
end
plot(x,M,'bo','MarkerFacecolor','b')
plot(xi,ML,'ro','MarkerFacecolor','r')

xi2=Eqd1dExact(x(1:(end-1)),Monitor);
ML2=interp1(x,M,xi2);
for ii = 1:4
  plot([xi2(ii) xi2(ii)],[ML2(ii) 0],'g-')
end
plot(xi2,ML2,'go','MarkerFacecolor','g')

title(['Equidistribution over a piecewise linear Monitor Function '...
    '\newline with 3 intervals (green) and 4 intervals (red).'])
ylabel('M')
xlabel('x')