function vis_permcut(ps)
global file_id
for js = 1:length(ps.tag)
    pd.dpi = 1;
    file_id = fopen(fullfile(ps.param.dir_source,[ps.tag{js},'.permcut.svg']),'w');
    fp = @(x) fprintf(file_id,[x,'\n']);
    ns = @(x) num2str(pd.dpi*x);
    
    sz = [180,45];
    maxsigma = 3;
    
    
    sigma = std(ps.z(~isnan(ps.z(:,js)),js));
    
    yp = ps.z(:,js)/sigma;
    xp = ps.shift;
    xpc = ps.shift(ps.icut{js});
    yp(abs(yp)>maxsigma)=maxsigma*sign(yp(abs(yp)>maxsigma));
    pd.size(1) = 800;
    putx = @(x) pd.size(1)*(5-mod(x,5))/5;
    pd.size(2) = 200;
    puty = @(x,y) pd.size(2)-(floor(x/5)*pd.size(2)/2 + pd.size(2)/4 + y/maxsigma*pd.size(2)/4);
    sy = pd.size(2)/2;
    fp(['<svg xmlns="http://www.w3.org/2000/svg" ' ...
        'width="',ns(pd.size(1)),'mm" height="',ns(pd.size(2)),'mm" ',...
        'viewBox="0 0 ',ns(pd.size(1)),' ',ns(pd.size(2)),'">']);
    svgo_rect([0,pd.size(1)],[0,pd.size(2)],'#ffffff');
    svgo_line([putx(0),putx(10)],[puty(1,0),puty(1,0)]);
    svgo_line([putx(0),putx(10)],[puty(9,0),puty(9,0)]);
    for i = 1:length(yp)
        svgo_line(putx(xp(i)),[puty(xp(i),0),puty(xp(i),yp(i))],'#000000',0.5);
    end
    for i = 1:length(xpc)
        svgo_line(putx(xpc(i)),[puty(xpc(i),-2),puty(xpc(i),2)],'#00ffff',0.5);
    end
    fp('</svg>');
    fclose(file_id);
end