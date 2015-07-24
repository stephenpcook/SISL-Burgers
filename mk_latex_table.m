function mystr = mk_latex_table(data,columnData,rowData,the_title,format,bold_lim,grey_lim)
% Makes a latex table of data where
% data(ii,jj) corresponds to columnData(ii) and rowData(jj)
% (so backwards from e.g. matrix indexing).
mystr = sprintf(['\\begin{table}[hbtp]\n \\caption{',the_title,'}']);
II = length(columnData);
JJ = length(rowData);
mystr = [mystr,sprintf('\n\\begin{tabular}{c||*{%d}{c|}} \n',II)];
mystr = [mystr,sprintf('& \\multicolumn{%d}{|c|}{N} \\\\ \\cline{2-%d} \n',II,II+1)];
mystr = [mystr,sprintf('$1/\\Delta t$')];
for ii = 1:II 
    mystr = [mystr,sprintf('& %d', columnData(ii))];
end % for ii
mystr = [mystr,sprintf('\\\\ \\hline \\hline \n')];
for jj = 1:JJ
  mystr = [mystr,sprintf('%d',rowData(jj))];
  for ii = 1:II 
    if data(ii,jj)>bold_lim
      mystr = [mystr,sprintf(['& {\\bf ',format,'}'], data(ii,jj))];
    elseif data(ii,jj)<grey_lim
      mystr = [mystr,sprintf(['& {\\color{gray} ',format,'}'], data(ii,jj))];
    else
      mystr = [mystr,sprintf(['& ',format], data(ii,jj))];
    end % if data(ii,jj)
  end % for ii
  mystr = [mystr,sprintf('\\\\ \n')];
end % for jj
mystr = [mystr,sprintf('\\hline \n\\end{tabular} \n\\end{table}')];
end % function mk_latex_table