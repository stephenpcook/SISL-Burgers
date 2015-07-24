%close all
%progs = {'burgersSLMM','burg2','burgersTom'};
%progs = {'burgersSLMM_MM'};
progs = {'burg2'}
%filename_end = '.mat';
filename_end = '_HRx.mat';

for jj = 1:length(progs)
figure
prog = progs{jj};
load(['out_',prog,filename_end])
loglog(bigN,(2/c*bigEps2(:,1)+1.5*c/bigTN(1)).*bigN','k-')
hold on
xlabel 'N_x'
ylabel 'K'
for ii = 2:length(bigTN)-1
loglog(bigN,(2/c*bigEps2(:,ii)+1.5*c/bigTN(ii)).*bigN','b-')
end % for ii
loglog(bigN,(2/c*bigEps2(:,end)+1.5*c/bigTN(end)).*bigN','r-')
hold off
title([prog,' : N*(bigEps2/0.5v + vDt) for different tN'])
xlabel 'N_x'
ylabel 'K'
end % for jj

for jj=1:length(progs)
figure
prog = progs{jj};
load(['out_',prog,filename_end])
loglog(bigN,bigEps2(:,1),'k-')
hold on
for ii=2:length(bigTN)-1
loglog(bigN,bigEps2(:,ii),'b-')
end % for ii
loglog(bigN,bigEps2(:,end),'r-')
title([prog,' : bigEps2 for different tN'])
xlabel 'N'
ylabel 'bigEps2'
end % for jj

for jj = 1:length(progs)
figure
prog = progs{jj};
load(['out_',prog,filename_end])
plot(bigN, bigC)
for ii = 1:length(bigTN)
  myLegend{ii} = ['N_t = ',num2str(bigTN(ii))];
end % for ii
legend(myLegend)
end % for jj