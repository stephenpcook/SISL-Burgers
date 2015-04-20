bigN = [20,40,80,160,320];
bigTN = [20,40,80,160,320]';

c = 1;
alpha = 0.5;
epsilon=0.0001;
grad = -alpha^2/(2*epsilon);

bigC = zeros(length(bigN),length(bigTN));
bigM = zeros(length(bigN),length(bigTN));
bigM2 = zeros(length(bigN),length(bigTN));


for ii = 1:length(bigN)
  for jj = 1:length(bigTN)
    N = bigN(ii);
    tN = bigTN(jj);
    %[U,X] = burgersSLMM(N,tN);
    [U,X] = burg2(N,tN);
    %[U,X] = burgersTom(N,tN);
    [m,cDt,w] = get_m_x(U,X,c,alpha);
    bigC(ii,jj) = cDt/1.5;
    bigM(ii,jj) = m;
    bigW(ii,jj) = w;
    %bigM2(ii,jj) = -m2;
  end % for jj
end % for ii

bigEps1 = -0.5*alpha^2./bigM;
bigEps2 = alpha*bigW/(4*atanh(0.95));

display(bigN)
display(bigTN)
display(bigC)

display(grad)
display(bigEps1)
display(bigEps2)

mystr = sprintf('\\begin{table} \n \\caption{Estimate of C}');
II = length(bigN);
JJ = length(bigTN);
mystr = [mystr,sprintf('\n\\begin{tabular}{c||*{%d}{c|}} \n',II)];
mystr = [mystr,sprintf('& \\multicolumn{%d}{|c|}{N} \\\\ \\cline{2-%d} \n',II,II+1)];
mystr = [mystr,sprintf('$1/\\Delta t$')];
for ii = 1:II 
    mystr = [mystr,sprintf('& %d', bigN(ii))];
end % for ii
mystr = [mystr,sprintf('\\\\ \\hline \\hline \n')];
for jj = 1:JJ
  mystr = [mystr,sprintf('%d',bigTN(jj))];
  for ii = 1:II 
    if bigC(ii)>1
      mystr = [mystr,sprintf('& {\\bf %.3f }', bigC(ii,jj))];
    else
      mystr = [mystr,sprintf('& %.3f', bigC(ii,jj))];
    end % if bigC
  end % for ii
  mystr = [mystr,sprintf('\\\\ \n')];
end % for jj
mystr = [mystr,sprintf('\\hline \n\\end{tabular} \n\\end{table}')];

mystr2 = sprintf('\\begin{table} \n \\caption{Estimate of $\\eps$ from gradient at $x=c$.}');
%II = length(bigN);
%JJ = length(bigTN);
mystr2 = [mystr2,sprintf('\n\\begin{tabular}{c||*{%d}{c|}} \n',II)];
mystr2 = [mystr2,sprintf('& \\multicolumn{%d}{|c|}{N} \\\\ \\cline{2-%d} \n',II,II+1)];
mystr2 = [mystr2,sprintf('$1/\\Delta t$')];
for ii = 1:II 
    mystr2 = [mystr2,sprintf('& %d', bigN(ii))];
end % for ii
mystr2 = [mystr2,sprintf('\\\\ \\hline \\hline \n')];
for jj = 1:JJ
  mystr2 = [mystr2,sprintf('%d',bigTN(jj))];
  for ii = 1:II 
    mystr2 = [mystr2,sprintf('& %.5f', -0.5*alpha^2/bigM(ii,jj))];
  end % for ii
  mystr2 = [mystr2,sprintf('\\\\ \n')];
end % for jj
mystr2 = [mystr2,sprintf('\\hline \n\\end{tabular}\n\\end{table}')];

mystr3 = sprintf('\\begin{table} \n \\caption{Estimate of $\\eps$ from width of front.}');
II = length(bigN);
JJ = length(bigTN);
mystr3 = [mystr3,sprintf('\n\\begin{tabular}{c||*{%d}{c|}} \n',II)];
mystr3 = [mystr3,sprintf('& \\multicolumn{%d}{|c|}{N} \\\\ \\cline{2-%d} \n',II,II+1)];
mystr3 = [mystr3,sprintf('$1/\\Delta t$')];
for ii = 1:II 
    mystr3 = [mystr3,sprintf('& %d', bigN(ii))];
end % for ii
mystr3 = [mystr3,sprintf('\\\\ \\hline \\hline \n')];
for jj = 1:JJ
  mystr3 = [mystr3,sprintf('%d',bigTN(jj))];
  for ii = 1:II 
    mystr3 = [mystr3,sprintf('& %.5f ', bigW(ii,jj)*(alpha/(4*1.832)))];
  end % for ii
  mystr3 = [mystr3,sprintf('\\\\ \n')];
end % for jj
mystr3 = [mystr3,sprintf('\\hline \n\\end{tabular} \n\\end{table}')];
