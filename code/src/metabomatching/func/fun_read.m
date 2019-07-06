function pr = fun_read(fn,format,hl)
if nargin<3
    hl=0;
end
fi = fopen(fn);
pr = textscan(fi,format,'delimiter','\t','HeaderLines',hl);
fclose(fi);
