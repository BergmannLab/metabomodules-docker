function svgo_text(x,y,txt,horiz,fsize,fweight,fcolor,fstyle)
% SVGO_TEXT  Write SVG object: text
global fontsize file_id
if nargin<8; fstyle='normal';  end
if nargin<7; fcolor='black';   end
if nargin<6; fweight='normal'; end
if nargin<5; fsize=fontsize;   end
if nargin<4; horiz='middle';   end

str = 'text ';
if round(x)==x && round(y)==y
    str = [str,sprintf('x="%03d" y="%03d"',x,y)];
else
    str = [str,sprintf('x="%06.2f" y="%06.2f"',x,y)];
end
str = [str,sprintf(' font-size="%02d"',fsize)];
str = [str,' font-family="Open Sans"'];
if ~strcmp(fweight,'normal'); str = [str,sprintf(' font-weight="%s"',fweight)]; end
if ~strcmp(fstyle, 'normal'); str = [str,sprintf(' font-style="%s"' ,fstyle )]; end
if ~strcmp(fcolor, 'black' ); str = [str,sprintf(' fill="%s"'       ,fcolor )]; end
if ~strcmp(horiz,  's'     ); str = [str,sprintf(' text-anchor="%s"',horiz  )]; end
fprintf(file_id,'<%s>%s</text>\n',str,txt);
