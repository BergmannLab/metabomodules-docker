function y = cas2num(x)
A=regexp(x,'-','split');
B=cell2mat(cellfun(@str2double,A,'uniformoutput',false));
C=B(:,1)*1e3;
D=B(:,2)*1e1;
E=B(:,3);
y=C+D+E;