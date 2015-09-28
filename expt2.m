% TODO - Make this a function that takes the in_filename as an input.

%close all
%progs = {'burgersSLMM','burg2','burgersTom'};
%progs = {'burgersSLMM_MM'};
program_name = 'burg2';
in_filename = 'test.mat';
load(in_filename)


for ii = 1:length(bigNT)
  myLegend{ii} = ['N_t = ',num2str(bigNT(ii))];
end % for ii


figure
loglog(bigNX,(2/c*bigEps1(:,1)+1.5*c/bigNT(1)).*bigNX','k-')
hold on
xlabel 'N_x'
ylabel 'K'
for ii = 2:length(bigNT)-1
loglog(bigNX,(2/c*bigEps1(:,ii)+1.5*c/bigNT(ii)).*bigNX','b-')
end % for ii
loglog(bigNX,(2/c*bigEps1(:,end)+1.5*c/bigNT(end)).*bigNX','r-')
hold off
title([program_name,' : N*(bigEps1/0.5v + vDt) for different tN'])
xlabel 'N_x'
ylabel 'K'
legend(myLegend)


figure
loglog(bigNX,bigEps1(:,1),'k-')
hold on
for ii=2:length(bigNT)-1
loglog(bigNX,bigEps1(:,ii),'b-')
end % for ii
loglog(bigNX,bigEps1(:,end),'r-')
title([program_name,' : bigEps1 for different tN'])
xlabel 'N_x'
ylabel 'bigEps1'
legend(myLegend)

figure
plot(bigNX, bigC)
xlabel 'N_x'
ylabel 'bigC'
legend(myLegend)