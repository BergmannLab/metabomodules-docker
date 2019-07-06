function tas=fun_readtable(fn,delimiter)

if nargin<2
    delimiter='\t';
end

fi = fopen(fn);
qq = textscan(fi,'%s',1,'delimiter','?');
qq = qq{1}{1};
fclose(fi);
nc = length(regexp(qq,'\s','split'));
fi = fopen(fn);
tas.lb = textscan(fi,repmat('%s',1,nc),1,'delimiter',delimiter);
tas.lb = [tas.lb{:}];
tas.lb = lower(tas.lb);
tas.pr = textscan(fi,repmat('%f',1,nc)  ,'delimiter',delimiter);
fclose(fi);
