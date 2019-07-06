function ps=vis_metabomatching(dir_source)
% VIS_METABOMATCHING  Create SVG images for metabomatching results
% dir_source=ps.param.dir_source;
ts.howto = false;
%% ##### COLORS #####
colhex.blue.darkBrewer   = '#1F78B4';
colhex.orange.darkBrewer = '#FF7F00';
colhex.green.darkBrewer  = '#33A02C';
colhex.green.liteBrewer  = '#B2DF8A';
colhex.green.darkBrewer  = '#006D2C';
colhex.green.middBrewer  = '#41AB5D';
colhex.green.liteBrewer  = '#C7E9C0';
colhex.purple.darkBrewer = '#54278F';
colhex.purple.middBrewer = '#807DBA';
colhex.purple.liteBrewer = '#DADAEB';
colhex.gray.liteBrewer   = '#D9D9D9';
colhex.gray.middBrewer   = '#969696';
%% ##### UNICODE KEYS #####
unicodeRep={'alpha','&#945;';'beta','&#946;';'gamma','&#947;';'delta','&#948;';'epsilon','&#949;';'omega','&#969;'};
%% ##### FILES ######
fn_parameters  = fullfile(dir_source,'parameters.out.tsv');
fn_metabolites = fullfile(dir_source,'metdb.mat');
fn_description = fullfile(dir_source,'description.tsv');
fn_control     = fullfile(dir_source,'cascontrol.tsv');
%% ##### GET PARAMETERS #####
ps = fun_load_parameters(dir_source,fn_parameters);
ts.scoreadj = ps.param.n_permutation>0;
% --- 2-compound mode implies many graphical changes
ts.is2c = ismember(ps.param.variant,{'2c','pm2c'});
% --- wider figure
if ismember('wide',fieldnames(ps.param))
    ts.wide = logical(ps.param.wide);
else
    ts.wide = false;
end
%% ##### GET PEAKS #####
if exist(fn_metabolites,'file');
    load(fn_metabolites);
else
    build_spectrum_database(dir_source);
end
%% ##### SVG FRIENDLY NAMES #####
rep_unicode={ ...
    'alpha-',  '&#945;-';'beta-' ,'&#946;-';...
    'gamma-',  '&#947;-';'delta-','&#948;-';...
    'epsilon-','&#949;-';'omega-','&#969;-'};
for k = 1:size(rep_unicode,1)
    metdb.name = strrep(metdb.name,rep_unicode{k,1},rep_unicode{k,2});
end
% ##### invalid CAS numbers #####
se = find(not(cellfun(@isempty,strfind(metdb.cas,'BMRB'))));
metdb.cas(se)={'0-0-0'};
metdb.casnum=cas2num(metdb.cas);
%% ##### GET PSEUDOSPECTRUM #####
ps = fun_import_pseudospectra(ps);
ps = fun_processps(ps);
%% ##### SET DEEP PARAMETERS #####
ps.param.deep.tag=strrep(strrep(ps.tag,'.','_'),'-','_');
pdeepdef={'ftt',true;'psa','left'};
for j = 1:size(pdeepdef,1)
    for i = 1:length(ps.tag)
        if ~ismember([pdeepdef{j,1},'_',ps.param.deep.tag{i}],fieldnames(ps.param));
            ps.param.deep.([pdeepdef{j,1},'_',ps.param.deep.tag{i}])=pdeepdef{j,2};
        else
            ps.param.deep.([pdeepdef{j,1},'_',ps.param.deep.tag{i}])=...
                ps.param.([pdeepdef{j,1},'_',ps.param.deep.tag{i}]);
        end
    end
end
switch ps.param.plot_type
    case 'p'
        ps.param.y_break_vl = 12;
        ps.param.y_majstep = 4;
        ps.param.y_minstep = 1;
        ps.param.y_break_pt = 3/4;
        ps.param.y_gmax = 256;
        ps.param.y_logset = [16,32:32:ps.param.y_gmax];
        ps.param.ps_ylabel = '&#x2212; log <tspan font-style="italic">p</tspan>';
    otherwise
        ps.param.y_break_vl = 8;
        ps.param.y_majstep = 2;
        ps.param.y_minstep = 1;
        ps.param.y_break_pt = 2/3;
        ps.param.y_gmax = 32;
        ps.param.y_logset = [8:8:ps.param.y_gmax];
        ps.param.ps_ylabel = 'Z-score';
end

if exist([dir_source,'op.csv'],'file')
    ps.op = csvread([dir_source,'op.csv']);
end
%% ##### CLEAN TAGS #####
for j = 1:length(ps.tag);
    tc = strrep(ps.tag{j},'.','_');
    tc = strrep(tc,'_neg','');
    tc = strrep(tc,'_pos','');
    ps.tac{j,1} = tc;
end
%% ##### GET DESCRIPTION
if exist(fn_description,'file')
    % description file can have 2 or 3 columns
    fi = fopen(fn_description);
    pr = textscan(fi,'%s',1,'delimiter','?');
    fclose(fi);
    xx = regexp(pr{1},'\t','split');
    nc = length(xx{1});
    if nc == 3
        fi = fopen(fn_description);
        pr = textscan(fi,'%s%s%s','delimiter','\t');
        fclose(fi);
        for i = 1:length(ps.tag)
            ix = find(strcmp(ps.tac{i},strrep(pr{1},'.','_')));
            if isempty(ix)
                ps.description{i,1}=['(',ps.tag{i},')'];
            else
                ps.description{i,1}=...
                    ['of ',pr{2}{ix},' in <tspan font-style="italic">',pr{3}{ix},'</tspan>'];
            end
        end
    else
        fi = fopen(fn_description);
        pr = textscan(fi,'%s%s','delimiter','\t');
        fclose(fi);
        for i = 1:length(ps.tag)
            ix = find(strcmp(ps.tac{i},strrep(pr{1},'.','_')));
            if isempty(ix)
                ps.description{i,1}='';
            else
                ps.description{i,1}=['of ',pr{2}{ix}];
            end
        end
    end
else
    for i = 1:length(ps.tag)
        if isfield(ps.param,['description_',ps.tac{i}])
            ps.description{i,1}=ps.param.(['description_',ps.tac{i}]);
        else
            ps.description{i,1}=['(',ps.tag{i},')'];
        end
    end
end
%% ##### GET CONTROL METABOLITES
for i = 1:length(ps.tag)
    ps.cas_control{i}={};
    ps.cas_control_num{i}=[];
end
if exist(fn_control,'file');
    fi = fopen(fn_control);
    pr = textscan(fi,'%s%s','delimiter','\t');
    fclose(fi);
    for i = 1:length(pr{1})
        field = pr{1}{i};
        value = pr{2}{i};
        field = strrep(field,'.','_');
        se = find(strcmp(field,ps.tac));
        if ~isempty(se)
            cas=value;
            casnum=cas2num(cas);
            for j = 1:length(se)
                ps.cas_control{se(j)}=[ps.cas_control{se(j)};cas];
                ps.cas_control_num{se(j)}=...
                    [ps.cas_control_num{se(j)};casnum];
            end
        end
    end
else
    for i = 1:length(ps.tag)
        aa = fieldnames(ps.param);
        se = find(strncmp(['cas_control_',ps.tac{i}],aa,12+length(ps.tac{i})));
        if not(isempty(se))
            for j = 1:length(se)
                ps.cas_control{i}{j}=ps.param.(aa{se(j)});
            end
        end
    end
end
%% ##### GET SCORES #####
clear matches;
if isfield(ps.param,'multi')
    fi = fopen(fullfile(dir_source,strrep(ps.param.multi,'pseudospectrum','score')));
    if ts.is2c
        fmt_met_lb = '%s%s%s%s';
        fmt_met_pr = '%s%f%s%f';
        nVar = 2;
    else
        fmt_met_lb = '%s%s';
        fmt_met_pr = '%s%f';
        nVar = 1;
    end
    lb = textscan(fi,[fmt_met_lb,repmat('%s',1,length(ps.tag))],1,'delimiter','\t');
    pr = textscan(fi,[fmt_met_pr,repmat('%f',1,length(ps.tag))],  'delimiter','\t');
    for J = 1:nVar
        matches.cas(:,J)=pr{(J-1)*2+1};
        matches.sid(:,J)=pr{(J-1)*2+2};
        [a,b] = ismember(matches.sid(:,J),metdb.sid);
        matches.name(a,J) = metdb.name(b(a));
        matches.casnum(a,J) = metdb.casnum(b(a));
    end
    for i = 1:length(ps.tag)
        xx = regexp(lb{2*nVar+i},'/','split');
        tag = xx{1}{2};
        ix = strcmpi(ps.tag,tag);
        matches.score(:,ix)=pr{2*nVar+i};
    end
    fclose(fi);
    if ts.scoreadj
        fi = fopen(fullfile(dir_source,strrep(ps.param.multi,'pseudospectrum','scoreadj')));
        if ts.is2c
            fmt_met_lb = '%s%s%s%s';
            fmt_met_pr = '%s%f%s%f';
            nVar = 2;
        else
            fmt_met_lb = '%s%s';
            fmt_met_pr = '%s%f';
            nVar = 1;
        end
        lb = textscan(fi,[fmt_met_lb,repmat('%s',1,length(ps.tag))],1,'delimiter','\t');
        pr = textscan(fi,[fmt_met_pr,repmat('%f',1,length(ps.tag))],  'delimiter','\t');
        for J = 1:nVar
            matches.cas(:,J)=pr{(J-1)*2+1};
            matches.sid(:,J)=pr{(J-1)*2+2};
            [a,b] = ismember(matches.sid(:,J),metdb.sid);
            matches.name(a,J) = metdb.name(b(a));
            matches.casnum(a,J) = metdb.casnum(b(a));
        end
        for i = 1:length(ps.tag)
            xx = regexp(lb{2*nVar+i},'/','split');
            tag = xx{1}{2};
            ix = strcmpi(ps.tag,tag);
            matches.scoreadj(:,ix)=pr{2*nVar+i};
        end
        fclose(fi);
    end
else
    for i = 1:length(ps.tag)
        fi = fopen(fullfile(dir_source,[ps.tag{i},'.score.tsv']));
        if ts.is2c
            pr = textscan(fi,'%s%f%s%f%f','delimiter','\t','headerlines',1);
            fclose(fi);
            if i == 1
                matches.cas(:,1)=pr{1};
                matches.sid(:,1)=pr{2};
                [a,b] = ismember(matches.sid(:,1),metdb.sid);
                matches.name(a,1) = metdb.name(b(a));
                matches.casnum(a,1) = metdb.casnum(b(a));
                matches.cas(:,2)=pr{3};
                matches.sid(:,2)=pr{4};
                [a,b] = ismember(matches.sid(:,2),metdb.sid);
                matches.name(a,2) = metdb.name(b(a));
                matches.casnum(a,2) = metdb.casnum(b(a));
            end
            matches.score(:,i)=pr{5};
            nVar=2;
        else
            pr = textscan(fi,'%s%f%f','delimiter','\t','headerlines',1);
            fclose(fi);
            if i == 1
                matches.cas=pr{1};
                matches.sid=pr{2};
                [a,b] = ismember(matches.sid,metdb.sid);
                matches.name = metdb.name(b(a));
                matches.casnum = metdb.casnum(b(a));
            end
            matches.score(:,i)=pr{3};
            nVar=1;
        end
        if ts.scoreadj
            fi = fopen(fullfile(dir_source,[ps.tag{i},'.scoreadj.tsv']));
            pr = textscan(fi,'%s%f%f','delimiter','\t','headerlines',1);
            fclose(fi);
            if i == 1
                matches.cas=pr{1};
                matches.sid=pr{2};
                [a,b] = ismember(matches.sid,metdb.sid);
                matches.name = metdb.name(b(a));
                matches.casnum = metdb.casnum(b(a));
            end
            matches.scoreadj(:,i)=pr{3};
            nVar=1;
        end
    end
end
%% ##### FORMAT PARAMETERS FOR SETTINGS BOX #####
%
% ----- figure title
pd.fg_set={};
%       pstype
if strcmp(ps.param.pstype,'correlation')
    pd.fg_set = [pd.fg_set;{'type','Correlation'}];
elseif strcmp(ps.param.pstype,'z')
    pd.fg_set = [pd.fg_set;{'type','Z-Score'}];
elseif strcmp(ps.param.pstype,'isa')
    pd.fg_set = [pd.fg_set;{'type','ISA'}];
elseif strcmp(ps.param.pstype,'pca')
    pd.fg_set = [pd.fg_set;{'type','PCA'}];    
end
%       variant
opt_variant = {'2-compound','&#x00B1;','&#x00B1;-2-compound'};
opt_ix = find(strcmp(ps.param.variant,{'2c','pm','pm2c'}));
if ~isempty(opt_ix)
    pd.fg_set = [pd.fg_set;{'variant',opt_variant{opt_ix}}];
end
%       mode
if strcmp(ps.param.mode,'peak')
    pd.fg_set = [pd.fg_set;{'mode',['peak, &#x03B4; = ',num2str(ps.param.dsh,'%.3f')]}];
elseif strcmp(ps.param.mode,'multiplet')
    pd.fg_set = [pd.fg_set;{'mode',['multiplet, &#x03B3; = ',num2str(ps.param.dsh,'%.3f')]}];
end
%       scoring
if strcmp(ps.param.scoring,'chisq')
    pd.fg_set = [pd.fg_set;{'scoring','&#x03C7;<tspan dy="-3" font-size="7">2</tspan>'}];
elseif strcmp(ps.param.scoring,'z')
    pd.fg_set = [pd.fg_set;{'scoring','Z'}];
end
%       reference database
if strcmp(ps.param.reference,'existing')
    pd.fg_set = [pd.fg_set;{'database','existing'}];
else
    pd.fg_set = [pd.fg_set;{'database',upper(ps.param.reference)}];
end
%       decorrelation
if ps.param.decorr_lambda == 0
    pd.fg_set = [pd.fg_set;{'decorr','&#x03BB; = 0'}];
elseif ps.param.decorr_lambda < 0.1
    pd.fg_set = [pd.fg_set;{'decorr',['&#x03BB; = ',num2str(ps.param.decorr_lambda,'%.2e')]}];
elseif ps.param.decorr_lambda < 1.0
    pd.fg_set = [pd.fg_set;{'decorr',['&#x03BB; = ',num2str(ps.param.decorr_lambda,'%.2f')]}];
end
%       decorrelation requested, but no correlation file provided
if isfield(ps.param,'correlation_missing')
    if ps.param.correlation_missing==1
        pd.fg_set = [pd.fg_set;{'warning','correlation file'}];
        pd.fg_set = [pd.fg_set;{'','missing'}];
    end
end
%% ##### ##### LOOPING OVER PSEUDOS ##### #####
for jPseudo=1:length(ps.tag)
    param=ps.param;
    param.description = ps.description{jPseudo};
    param.cas_ctrl = ps.cas_control{jPseudo};
    param.tag = ps.tag{jPseudo};
    pseudo.x = ps.shift;
    switch param.plot_type
        case 'p'
            pseudo.beta = ps.beta(:,jPseudo);
            pseudo.c = 1+(pseudo.beta<0);
            pseudo.y = -log10(ps.p(:,jPseudo));
            % pseudo.x_sig = pseudo.y>=param.significant & abs(pseudo.beta)>1e-5;
        otherwise
            pseudo.z = ps.z(:,jPseudo);
            pseudo.c = 1+(pseudo.z<0);
            pseudo.y = abs(pseudo.z);
    end
    if ~isfield(param,'significant')
        param.significant = max(pseudo.y);
        if isnan(param.significant)
            break
        end
    end
    pseudo.x_sig = pseudo.y>=param.significant;
    
    match = matches;
    match.score = matches.score(:,jPseudo);
    if ts.scoreadj
        match.scoreadj=matches.scoreadj(:,jPseudo);
    end
    %
    % ----- pseudospectrum title
    pd.fg_pseudo_title = ['Pseudospectrum ',param.description];
    %       figure requires number of candidates to be a multiple of 4
    nshow = 4*round(param.n_show/4);
    %
    % ------+----------------------------------------------+
    %       | COMPUTE RANKS AND DEFINE METABOLITES TO SHOW |
    %       +----------------------------------------------+---------------
    % ----- compute rank
    match.score(isnan(match.score))=0;
    [~,A] = sort(match.score,'descend');
    [~,B] = sort(A,'ascend');
    match.rank = B;
    [a,b] = sort(match.rank);
    for f = fieldnames(match)'
        match.(f{1})=match.(f{1})(b,:);
    end
    % ----- check whether subsequent spectra belong to the same metabolite
    %       or metabolite pair; if so, group them
    A1=match.casnum(2:end,:);
    A2=match.casnum(1:(end-1),:);
    TT=all(A1==A2,2);
    if ts.is2c
        S1=A1(:,1)==A2(:,2);
        S2=A1(:,2)==A2(:,1);
        TT=TT|(S1&S2);
    end
    for f = fieldnames(match)'
        match.(f{1})(TT,:)=[];
    end
    % ----- assing ranks, define what to show
    match.rank=(1:length(match.rank))';
    match.show = 0*match.rank;
    match.ctrl = 0*match.rank;
    sl = match.rank<=nshow;
    match.show(sl)=match.rank(sl);
    se_ctrl=[];
    if ~isempty(ps.cas_control)
        if ts.is2c
            if iscell(param.cas_ctrl)
                jc=0;
                for ic = 1:length(param.cas_ctrl)
                    a = find(strcmp(match.cas(:,1),param.cas_ctrl{ic}));
                    if isempty(a)
                        fi=fopen(fullfile(dir_source,'info.log'),'a');
                        fprintf(fi,'Provided CASRN %s for control metabolite does not exist in reference database',param.cas_ctrl{ic});
                        fclose(fi);
                    end
                    b = find(strcmp(match.cas(:,2),param.cas_ctrl{ic}));
                    se = unique([a;b]);
                    [~,ix] = min(match.rank(se));
                    if isempty(ix)
                        se_ctrl=unique([se_ctrl;se(ix)]);
                    end
                end
                if length(param.cas_ctrl)>1
                    c = find(strcmp(match.cas(:,1),param.cas_ctrl{1})&strcmp(match.cas(:,2),param.cas_ctrl{2}));
                    d = find(strcmp(match.cas(:,2),param.cas_ctrl{1})&strcmp(match.cas(:,1),param.cas_ctrl{2}));
                    se = unique([c;d]);
                    [~,ix]=min(match.rank(se));
                    if ~isempty(ix)
                        se_ctrl=unique([se_ctrl;se(ix)]);
                    end
                end
            else
                a = find(strcmp(match.cas(:,1),param.cas_ctrl));
                b = find(strcmp(match.cas(:,2),param.cas_ctrl));
                se = unique([a;b]);
                [~,ix] = min(match.rank(se));
                if ~isempty(ix)
                    se_ctrl=se(ix);
                end
            end
            
        else
            if iscell(param.cas_ctrl)
                for i = 1:length(param.cas_ctrl)
                    ix = find(strcmp(match.cas,param.cas_ctrl{i}));
                    if ~isempty(ix)
                        se_ctrl=[se_ctrl;ix];
                    end
                end
            else
                se_ctrl = find(strcmp(match.cas,param.cas_ctrl));
            end
        end
    end
    for i=1:length(se_ctrl)
        ix = se_ctrl(i);
        match.ctrl(ix)=1;
        if match.show(ix)==0
            match.show(ix)=match.rank(ix);
        end
    end
    while sum(match.show>0)>nshow
        se = find(match.show>0 & match.ctrl==0);
        [~,ix] = max(match.rank(se));
        match.show(se(ix))=0;
    end
    se=find(match.show>0);
    [~,reo0]=sort(match.rank(se));
    match.show(se(reo0))=(1:nshow)';
    
    pd.xRange = [0,10];
    pd.xDiff = pd.xRange(2)-pd.xRange(1);
    
    A={colhex.blue.darkBrewer,colhex.orange.darkBrewer};
    pseudo.c = A(pseudo.c);
    [~,reo] = sort(pseudo.x);
    for f={'x','y','c','x_sig'}
        pseudo.(f{1})=pseudo.(f{1})(reo);
    end
    
    for i = 1:nshow
        ix = find(match.show==i);
        spectrum.score(i) = match.score(ix);
        if ts.scoreadj
            spectrum.scoreadj(i) = match.scoreadj(ix);
        end
        for j = 1:nVar
            jx = find(metdb.sid==match.sid(ix,j));
            spectrum.met_name{i,j} = metdb.name{jx};
            spectrum.met_cas{i,j} = match.cas{ix,j};
            if strcmp(param.mode,'peak')
                spectrum.met_pos1{i,j} = metdb.pos{jx}-param.dsh;
                spectrum.met_pos2{i,j} = metdb.pos{jx}+param.dsh;
            elseif strcmp(param.mode,'multiplet');
                spectrum.met_pos1{i,j} = metdb.pos{jx}(:,1)-param.dsh;
                spectrum.met_pos2{i,j} = metdb.pos{jx}(:,2)+param.dsh;
            end
            spectrum.met_h{i,j} = metdb.hnh{jx}/max(metdb.hnh{jx});
        end
        spectrum.rank(i) = match.rank(ix);
        spectrum.show(i) = match.show(ix);
        spectrum.ctrl(i) = match.ctrl(ix);
    end
    
    % significant peak clusters
    
    sh_cl=[];
    sh_sl=[];
    sh_label=[];
    sh_p=[];
    
    if any(pseudo.x_sig)
        sh=pseudo.x(pseudo.x_sig);
        sh_p = pseudo.y(pseudo.x_sig);
        
        sh_cl=1;
        ic=1;
        if length(sh)>1
            for i = 2:length(sh)
                if sh(i)-sh(i-1)>0.03;
                    ic=ic+1;
                end
                sh_cl(i)=ic;
                
            end
        end
        for i = 1:max(sh_cl)
            sh_sl = sh(sh_cl==i);
            p_sl=sh_p(sh_cl==i);
            [a,ix]=min(p_sl);
            sh_label(i)=sh_sl(ix);
        end
    end
    % scale p-values
    %     if all(pseudo.y<1)
    %         pseudo.y = -log10(pseudo.y);
    %     end
    pseudo.y(end+1) = param.significant; % adding p_sig here to get its val in new axis
    pseudo.y(pseudo.y>param.y_gmax)=param.y_gmax;
    sl1=pseudo.y< param.y_break_vl;
    sl2=pseudo.y>=param.y_break_vl;
    pseudo.y(sl1)=param.y_break_pt*pseudo.y(sl1)/param.y_break_vl;
    pseudo.y(sl2)=param.y_break_pt+(1-param.y_break_pt)*...
        log2(1+pseudo.y(sl2)-param.y_break_vl)/log2(1+param.y_gmax);
    pseudo.ySig=pseudo.y(end);
    pseudo.y(end)=[];
    majset = (0:param.y_majstep:param.y_break_vl);
    minset = setdiff((0:param.y_minstep:param.y_break_vl),majset);
    pseudo.yax_maj = ...
        1-[param.y_break_pt*majset/param.y_break_vl,1];
    pseudo.yax_maj_vl = [majset,param.y_gmax];
    pseudo.yax_min = [...
        param.y_break_pt*minset/param.y_break_vl,...
        param.y_break_pt+(1-param.y_break_pt)*log2(1+param.y_logset)/log2(1+param.y_gmax)];
    pd.dpi = 1;
    
    pd.d0 = 4;
    pd.d1 = 4;
    pd.d_tmin = 3;
    pd.d_tmaj = 5;
    pd.d_text_hght = 8;
    pd.d_text_wdth = 5.5;
    pd.d_text_midd = 3;
    pd.d_text_spce = 4;
    pd.d_text_line=pd.d_text_hght+pd.d_text_spce;
    pd.d_spctr_sp = 5;
    pd.dl_grid = 0.5;
    pd.c_grid = colhex.gray.middBrewer;
    pd.c_mark = 'white';
    if isfield(ps.param,'mark')
        switch ps.param.mark
            case 'gray'
                pd.c_mark = colhex.gray.middBrewer;
        end
    end
    pd.dl_pseudo = 0.75;
    if ts.wide
        pd.size = [800,700];
        pd.n_maxchar = 25;
        pd.x1 = 520;
        pd.x2 = pd.size(1)-pd.x1;
        pd.d_siglab = 0.32; % minimum distance between labels of significant features (in ppm)
    else
        pd.size = [600,700];
        pd.n_maxchar = 15;
        pd.x1 = 440;
        pd.x2 = pd.size(1)-pd.x1;
        pd.d_siglab = 0.40; % minimum distance between labels of significant features (in ppm)
    end
    param.yPseudo = 100;
    pd.yyy = [40,100,30,100,100];
    pd.x_yaxis = 32;
    pd.yp = 0;
    
    pd.dFromTitle = 3*pd.d1;
    pd.d_pseudo_title = 2*pd.d1+pd.d_text_hght;
    pd.d_pseudo_boxxx = 2*pd.d1;
    pd.dFromSpectrumTitle = 2*pd.d1;
    pd.dFromSpectrum = 2*pd.d1;
    
    pd.o2=pd.x2+pd.d0;
    pd.s2=pd.x2-2*pd.d0;
    pd.s1=pd.x1-2*pd.d0;
    
    if ismember(['ftt_',ps.param.deep.tag{jPseudo}],fieldnames(ps.param.deep))
        ts.show_title=ps.param.deep.(['ftt_',ps.param.deep.tag{jPseudo}]);
    else
        ts.show_title = true;
    end
    if ts.show_title
        pd.o_fig_title = [pd.d0,pd.d0];
        pd.s_fig_title = [pd.x2-2*pd.d0,param.yPseudo];
        pd.yp=pd.d0;
    else
        pd.yp = pd.d0;
    end
    
    pd.yp = pd.d0;
    pd.o_pseudo_title = [pd.o2,                     pd.yp];
    pd.s_pseudo_title = [pd.s1,                     pd.d_text_hght+pd.d_text_spce];
    pd.o_pseudo_yaxis = [-pd.x_yaxis,               0];
    pd.s_pseudo_yaxis = [ pd.x_yaxis,               param.yPseudo];
    pd.yp = pd.o_pseudo_title(2)+pd.s_pseudo_title(2)+pd.d_pseudo_title;
    
    pd.o_pseudo_boxxx = [pd.o2+pd.s_pseudo_yaxis(1),pd.yp];
    pd.s_pseudo_boxxx = [pd.s1-pd.s_pseudo_yaxis(1),param.yPseudo];
    pd.o_pseudo_xaxis = [0,                         pd.s_pseudo_boxxx(2)];
    pd.s_pseudo_xaxis = [pd.s_pseudo_boxxx(1),      pd.yyy(3)];
    pd.yp = pd.o_pseudo_boxxx(2)+pd.s_pseudo_boxxx(2)+pd.s_pseudo_xaxis(2)+pd.d_pseudo_boxxx;
    
    pd.o_spectr_title = [pd.d0,                     pd.yp];
    pd.s_spectr_title = [pd.s2,                     pd.d_text_hght];
    pd.yp = pd.o_spectr_title(2)+pd.s_spectr_title(2)+pd.dFromSpectrumTitle;
    
    pd.o_spectr_subtt = [pd.o_pseudo_boxxx(1),      pd.yp];
    pd.s_spectr_subtt = [pd.s_pseudo_boxxx(1),      pd.d_text_hght];
    pd.o_spectr_headr = [pd.o_pseudo_yaxis(1),      0];
    pd.s_spectr_headr = [pd.s_pseudo_yaxis(1),      pd.d_text_hght];
    pd.yp = pd.o_spectr_subtt(2)+pd.s_spectr_subtt(2);
    
    pd.o_spectr_boxxx = [pd.o_pseudo_boxxx(1),      pd.yp+pd.d1];
    pd.y_spectr_lines = pd.d_spctr_sp/2+pd.d_text_line*(1:nshow)+(ceil((1:nshow)/4)-1)*pd.d_spctr_sp;
    pd.y_spectr_boxxx = pd.y_spectr_lines(end)+pd.d_spctr_sp/2;
    pd.s_spectr_boxxx = [pd.s_pseudo_boxxx(1),      pd.y_spectr_boxxx];
    pd.o_spectr_label = [pd.o_pseudo_yaxis(1),      0];
    pd.s_spectr_label = [pd.s_pseudo_yaxis(1),      pd.y_spectr_boxxx];
    pd.yp = pd.o_spectr_boxxx(2)+pd.s_spectr_boxxx(2)+pd.dFromSpectrum;
    
    pd.o_spectr_xaxis = [0,                         pd.s_spectr_boxxx(2)];
    pd.s_spectr_xaxis = [pd.o_pseudo_boxxx(1),      pd.yyy(3)];
    
    pd.yp = pd.o_spectr_boxxx(2)+pd.s_spectr_boxxx(2)+pd.s_spectr_xaxis(2)+pd.d_pseudo_boxxx;
    
    pd.o_metabo_title = [pd.d0,                     pd.o_spectr_subtt(2)];
    pd.s_metabo_title = [pd.s2,                     pd.s_spectr_subtt(2)];
    pd.yp = pd.o_metabo_title(2)+pd.s_metabo_title(2);
    
    pd.o_metabo_boxxx = [pd.d0,                     pd.yp+pd.d1];
    pd.s_metabo_boxxx = [pd.s_metabo_title(1),      pd.s_spectr_boxxx(2)];
    
    pd.yp = pd.o_metabo_boxxx(2)+pd.s_metabo_boxxx(2)+pd.s_spectr_xaxis(2)+pd.d0;
    
    pd.size(2) = pd.yp;
    
    ns = @(x) num2str(pd.dpi*x);
    
    sh2x = @(ppm) pd.s_pseudo_boxxx(1)*pd.dpi*(1-(ppm-pd.xRange(1))/pd.xDiff);
    add0 = @(str) ['<tspan fill="white">0</tspan>',str];
    
    %
    global file_id fontsize
    if ts.howto
        file_id = fopen(fullfile(param.dir_source,'f1.svg'),'w');
    else
        file_id = fopen(fullfile(param.dir_source,[param.tag,'.svg']),'w');
    end
    fontsize = 10;
    fp = @(x) fprintf(file_id,[x,'\n']);
    ggo = @(x) fp(['<g transform="translate(',num2str(x(1),'%03d'),',',num2str(x(2),'%03d'),')">']);
    %
    
    sz=[180,180/pd.size(1)*pd.size(2)];
    fp(['<svg xmlns="http://www.w3.org/2000/svg" ' ...
        'width="',ns(pd.size(1)),'mm" height="',ns(pd.size(2)),'mm" ',...
        'viewBox="0 0 ',ns(pd.size(1)),' ',ns(pd.size(2)),'">']);
    svgo_rect([0,pd.size(1)],[0,pd.size(2)],'#ffffff');
    % ##### TITLE GROUP #####
    if ts.show_title
        ggo(pd.o_fig_title);
        if ts.howto
            svgo_rect([-2,pd.s_fig_title(1)*.86],[-2,pd.s_fig_title(2)*.61],'#ffffff','#e31a1c');
            svgo_text_c(pd.s_fig_title(1)*.86+2,mean([-2,pd.s_fig_title(2)*.61])+3,'A','howto','start');
        end
        svgo_text_c(0,pd.d_text_hght,'Metabomatching Settings','heading','s');
        for j = 1:size(pd.fg_set,1)
            svgo_text_c( 0,pd.o_pseudo_boxxx(2)-1+(j-1)*pd.d_text_line,pd.fg_set{j,1},'normal','s');
            svgo_text_c(52,pd.o_pseudo_boxxx(2)-1+(j-1)*pd.d_text_line,pd.fg_set{j,2},'normal','s');
        end
        fp('</g>');
    end
    % ##### PSEUDOSPECTRUM TITLE GROUP #####
    ggo(pd.o_pseudo_title);
    if ts.howto
        svgo_rect([-4,sh2x(5.35)],[-2,pd.s_pseudo_title(2)],'#ffffff','#e31a1c');
        svgo_text_c(sh2x(5.35)+2,mean([-2,pd.s_pseudo_title(2)])+3,'B','howto','start');
    end
    if ismember(['psa_',ps.param.deep.tag{jPseudo}],fieldnames(ps.param.deep))
        if strcmp(ps.param.deep.(['psa_',ps.param.deep.tag{jPseudo}]),'right')
            svgo_text_c(pd.s_pseudo_title(1),pd.d_text_hght,pd.fg_pseudo_title,'heading','end');
        else
            svgo_text_c(0,pd.d_text_hght,pd.fg_pseudo_title,'heading','start');
        end
    else
        svgo_text_c(0,pd.d_text_hght,pd.ssTitlePseudo,'heading','start');
    end
    fp('</g>');
    %
    %
    % ##### PSEUDOSPECTRUM GROUP #####
    ggo(pd.o_pseudo_boxxx);
    if ts.howto
        svgo_rect([sh2x(3.4),sh2x(0.9)],[1,-3]*pd.d_text_hght,'#ffffff','#e31a1c');
        svgo_text_c(sh2x(0.9)+2,mean([1,-3]*pd.d_text_hght)+3,'C','howto','start');
        svgo_rect([sh2x(9.8),sh2x(10.9)],[-6,pd.s_pseudo_boxxx(2)+6],'#ffffff','#e31a1c');
        svgo_text_c(sh2x(10.9),pd.s_pseudo_boxxx(2)+17,'D','howto','start');
    end
    %
    % ----- grid lines
    %       vertical
    A=[0,pd.s_pseudo_boxxx(1),sh2x(1:9)];
    for ia=1:length(A)
        svgo_line(A(ia),[0,pd.s_pseudo_boxxx(2)],pd.c_grid,pd.dl_grid);
    end
    %       horizontal
    A=pd.s_pseudo_boxxx(2)*[0,(1-param.y_break_pt),1];
    for ia=1:length(A)
        svgo_line([0,pd.s_pseudo_boxxx(1)],A(ia),pd.c_grid,pd.dl_grid);
    end
    svgo_line([0,pd.s_pseudo_boxxx(1)],(1-pseudo.ySig)*pd.s_pseudo_boxxx(2),colhex.green.liteBrewer,pd.dl_grid);
    %
    % ----- markers for significant peaks
    Q=find(pseudo.x_sig);
    for i = 1:length(Q)
        svgo_line(sh2x(pseudo.x(Q(i))),pd.d1*[-1,1]/2,pd.c_grid,pd.dl_grid);
    end
    %       labels (one per cluster)
    ypos=[];
    sh_leftmost=[];
    ypos(1)=1;
    for i=1:length(sh_label)
        if i>1
            for j=1:max(ypos)
                find(ypos==j,1,'last');
                sh_leftmost(j) = sh_label(find(ypos==j,1,'last'));
            end
            sl = abs(sh_label(i)-sh_leftmost)>pd.d_siglab;
            ix = find(sl,1,'first');
            if isempty(ix)
                ypos(i)=max(ypos)+1;
            else
                ypos(i) = ix;
            end
        end
        
        nnn = num2str(round(100*sh_label(i)),'%03d');
        
        svgo_text_c(sh2x(sh_label(i)),-pd.d1-(ypos(i)-1)*pd.d_text_hght,[nnn(1),...
            '<tspan dy="-1" dx="-.5" font-size="8">',nnn(2),'</tspan>',...
            '<tspan dx="-.5" font-size="8">',nnn(3),'</tspan>'],'legend');
    end
    % ----- pseudospectrum
    for i = 1:length(pseudo.x)
        svgo_line(sh2x(pseudo.x(i)),[1,1-pseudo.y(i)]*pd.s_pseudo_boxxx(2),pseudo.c{i},pd.dl_pseudo);
    end
    % ----- pseudospectrum x-axis
    ggo(pd.o_pseudo_xaxis);
    if ts.howto
        svgo_rect([sh2x(1),sh2x(0)+3.5],[17,32],'#ffffff','#e31a1c');
        svgo_text_c(sh2x(1)-2,pd.d1+1.9*pd.d_text_hght+2*pd.d_text_spce+1,'E','howto','end');
    end
    for i = 0:0.2:10
        if i<=pd.xRange(2) && i>=pd.xRange(1)
            if mod(i,1)<1e-4
                svgo_line(sh2x(i),pd.d_tmaj*[-1,1],pd.c_grid,pd.dl_grid);
                svgo_text_c(sh2x(i),pd.d1+pd.d_text_hght+pd.d_text_spce*0.75,num2str(i),'normal');
            end
        end
    end
    svgo_text_c(sh2x(mean(pd.xRange)),pd.d1+2*pd.d_text_hght+2*pd.d_text_spce,'chemical shift [ppm]','normal');
    pd.ssBetaLegend = [...
        '&#946;:',...
        '<tspan dy="-0.5" fill="',colhex.orange.darkBrewer,  '">&#9632;</tspan>',...
        '<tspan dy="+0.5" dx="-0.3">&#60;0&#60;</tspan>',...
        '<tspan dy="-0.5" dx="-0.3" fill="',colhex.blue.darkBrewer,  '">&#9632;</tspan>'];
    svgo_text_c(3+sh2x(min(pd.xRange)),pd.d1+1.9*pd.d_text_hght+2*pd.d_text_spce,pd.ssBetaLegend,'legend','end');
    fp('</g>');
    % ----- pseudospectrum y-axis
    ggo(pd.o_pseudo_yaxis);
    svgo_line(pd.s_pseudo_yaxis(1),[0,pd.s_pseudo_yaxis(2)],pd.c_grid,pd.dl_grid);
    fprintf(file_id,[...
        '<text x="%f" y="%f" ',...
        'font-size="',num2str(fontsize),...
        '" font-family="Open Sans',...
        '" text-anchor="middle',...
        '" transform="rotate(%d,%f,%f)',...
        '">%s</text>\n'],...
        pd.d_text_hght,pd.s_pseudo_yaxis(2)/2,270,pd.d_text_hght,pd.s_pseudo_yaxis(2)/2,param.ps_ylabel);
    for i = 1:length(pseudo.yax_maj)
        svgo_line(pd.s_pseudo_yaxis(1)-pd.d_tmaj*[1,-1],pd.s_pseudo_yaxis(2)*pseudo.yax_maj(i),pd.c_grid,pd.dl_grid);
        svgo_text_c(pd.s_pseudo_yaxis(1)-pd.d_tmaj,pd.s_pseudo_yaxis(2)*pseudo.yax_maj(i)+3,num2str(pseudo.yax_maj_vl(i)),'normal','end');
    end
    for i = 1:length(pseudo.yax_min)
        svgo_line(pd.s_pseudo_yaxis(1)-pd.d_tmin*[-1,0],pd.s_pseudo_yaxis(2)*(1-pseudo.yax_min(i)),pd.c_grid,pd.dl_grid);
    end
    fp('</g>');
    fp('</g>');
    % ##### SPECTRUM TITLE GROUP #####
    ggo(pd.o_spectr_title);
    svgo_text_c(0,pd.d_text_hght,'Candidate Metabolites','heading','start');
    fp('</g>');
    % ##### SPECTRUM SUBTITLE GROUP #####
    ggo(pd.o_spectr_subtt);
    svgo_text_c(0,pd.d_text_hght,'match set','normal','s');
    ggo(pd.o_spectr_headr);
    if ts.scoreadj
        svgo_text_c(0,pd.d_text_hght,'adj','normal','s');
    else
        svgo_text_c(0,pd.d_text_hght,'rank','normal','s');
    end
    fp('</g>');
    fp('</g>');
    % ##### SPECTRUM GROUP #####
    ggo(pd.o_spectr_boxxx);
    if ts.howto
        svgo_rect(sh2x([1.1,2.75]),pd.y_spectr_lines(2)-pd.d_text_line/2+pd.d_text_hght*6/9*[-1,1],'#ffffff','#e31a1c');
        svgo_text_c(sh2x(2.75)-2,pd.y_spectr_lines(2)-2,'G','howto','end');
        % svgo_rect([-pd.s_pseudo_yaxis(1),-16]-2,pd.y_spectr_lines(1)-pd.d_text_line/2+pd.d_text_hght*6/9*[-1,1],'#ffffff','#e31a1c');
        % svgo_text_c(-pd.s_pseudo_yaxis(1)-5,pd.y_spectr_lines(1)-2,'F','howto','end');
        
    end
    svgo_line(0,[0,pd.s_spectr_boxxx(2)],pd.c_grid,pd.dl_grid);
    svgo_line(pd.s_spectr_boxxx(1),[0,pd.s_spectr_boxxx(2)],pd.c_grid,pd.dl_grid);
    Z=nshow/4;
    for i=1:(Z-1)
        svgo_line([-pd.s_pseudo_yaxis(1),pd.s_spectr_boxxx(1)],mean(pd.y_spectr_lines([4*i,4*i+1]))-pd.d_text_line/2,pd.c_grid,pd.dl_grid);
    end
    svgo_line([-pd.s_pseudo_yaxis(1),pd.s_spectr_boxxx(1)],0,pd.c_grid,pd.dl_grid);
    svgo_line([-pd.s_pseudo_yaxis(1),pd.s_spectr_boxxx(1)],pd.s_spectr_boxxx(2),pd.c_grid,pd.dl_grid);
    for i = 1:9
        svgo_line(sh2x(i),[0,pd.s_spectr_boxxx(2)],pd.c_grid,pd.dl_grid);
    end
    % ----- markers for significant peaks
    for i=1:length(sh_label)
        nnn = num2str(round(100*sh_label(i)),'%03d');
        svgo_text_c(sh2x(sh_label(i)),-pd.d1-(ypos(i)-1)*pd.d_text_hght,[nnn(1),...
            '<tspan dy="-1" dx="-.6" font-size="8">',nnn(2),'</tspan>',...
            '<tspan dx="-.6" font-size="8">',nnn(3),'</tspan>'],'legend');
    end
    Q=find(pseudo.x_sig);
    for i = 1:length(Q)
        for j = 1:Z-1
            svgo_line(sh2x(pseudo.x(Q(i))),mean(pd.y_spectr_lines([4*j,4*j+1]))-pd.d_text_line/2+pd.d1*[-1,1]/2,pd.c_grid,pd.dl_grid);
        end
        svgo_line(sh2x(pseudo.x(Q(i))),pd.d1*[-1,1]/2,pd.c_grid,pd.dl_grid);
        svgo_line(sh2x(pseudo.x(Q(i))),pd.s_spectr_boxxx(2)+pd.d1*[-1,1]/2,pd.c_grid,pd.dl_grid);
    end
    % ----- spectra peak rectangles
    cols = {colhex.green.darkBrewer,colhex.green.middBrewer,colhex.green.liteBrewer;...
        colhex.purple.darkBrewer,colhex.purple.middBrewer,colhex.purple.liteBrewer};
    if isfield(spectrum,'metaboliteColor'); spectrum=rmfield(spectrum,'metaboliteColor'); end
    for i = 1:nshow
        for j = 1:nVar
            for k = 1:length(spectrum.met_h{i,j})
                spectrum.metaboliteColor{i,j}{k,1} = '#000000';
                if spectrum.met_h{i,j}(k)>0.8
                    spectrum.metaboliteColor{i,j}{k,1} = cols{j,1};
                elseif spectrum.met_h{i,j}(k)>0.2
                    spectrum.metaboliteColor{i,j}{k,1} = cols{j,2};
                else
                    spectrum.metaboliteColor{i,j}{k,1} = cols{j,3};
                end
            end
        end
        
    end
    for i = 1:nshow
        hh = [];
        p1 = [];
        p2 = [];
        cc = {};
        qq = [];
        for j = 1:nVar
            hh=[hh;spectrum.met_h{i,j}];
            p1=[p1;spectrum.met_pos1{i,j}];
            p2=[p2;spectrum.met_pos2{i,j}];
            cc=[cc;spectrum.metaboliteColor{i,j}];
            qq=[qq;(-1)^j*ones(length(spectrum.met_h{i,j}),1)];
        end
        if nVar==1; qq=0*qq; end
        [~,reo]=sort(hh,'ascend');
        for j = reo'
            svgo_rect(sh2x([p1(j),p2(j)]),...
                pd.y_spectr_lines(i)-pd.d_text_line/2+qq(j)*pd.d_text_hght*1/9+pd.d_text_hght*(5-abs(qq(j)))/9*[-1,1],cc{j});
        end
    end
    for i = 1:nshow
        hh = [];
        p1 = [];
        p2 = [];
        cc = {};
        for j = 1:nVar
            hh=[hh;spectrum.met_h{i,j}];
            p1=[p1;spectrum.met_pos1{i,j}];
            p2=[p2;spectrum.met_pos2{i,j}];
            cc=[cc;spectrum.metaboliteColor{i,j}];
        end
        [~,reo]=sort(hh,'ascend');
        for j = reo'
            seTag=find(...
                p1(j)-pseudo.x(Q)<0 & ...
                pseudo.x(Q)-p2(j)<0);
            if ~isempty(seTag)
                for k=1:length(seTag)
                    svgo_line(sh2x(pseudo.x(Q(seTag(k)))),pd.y_spectr_lines(i)-pd.d_text_line/2+pd.d_text_hght*1/3*[-1,1],pd.c_mark,pd.dl_pseudo);
                end
            end
        end
    end
    % ##### SPECTRUM Y-LABEL BOX #####
    ggo(pd.o_spectr_label);
    % ----- spectrum rank
    %       if permutation scores exist, print them instead of ranks
    if ts.scoreadj
        for i = 1:length(pd.y_spectr_lines)
            str = num2str(spectrum.scoreadj(i),'%.1f');
            svgo_text_c(0,pd.y_spectr_lines(i)-pd.d_text_midd,str,'normal','s');
        end
    else
        rankExponent=max(floor(log10(spectrum.rank)));
        for i = 1:length(pd.y_spectr_lines)
            A = floor(log10(spectrum.rank(i)));
            str = num2str(spectrum.rank(i));
            while A < rankExponent
                str=add0(str);
                A=A+1;
            end
            svgo_text_c(0,pd.y_spectr_lines(i)-pd.d_text_midd,str,'normal','s');
        end
    end
    fp('</g>');
    % ----- spectrum x-axis
    ggo(pd.o_spectr_xaxis);
    if ts.howto
        svgo_rect([sh2x(2.9),sh2x(0)+3.5],[17,32],'#ffffff','#e31a1c');
        svgo_text_c(sh2x(2.9)-2,pd.d1+1.9*pd.d_text_hght+2*pd.d_text_spce+1,'H','howto','end');
    end
    for i = 0:0.2:10
        if i<=pd.xRange(2) && i>=pd.xRange(1)
            if mod(i,1)<1e-4
                svgo_line(sh2x(i),pd.d_tmaj*[-1,1],pd.c_grid,pd.dl_grid);
                svgo_text_c(sh2x(i),pd.d1+pd.d_text_hght+pd.d_text_spce*0.75,num2str(i),'normal');
            else
                %svgo_line(sh2x(i),pd.s_pseudo_xaxis(2)-[0,pd.d1/2]);
                %svgo_line(sh2x(i),pd.d_tmin*[-1,0]);
            end
        end
    end
    svgo_text_c(sh2x(mean(pd.xRange)),pd.d1+2*pd.d_text_hght+2*pd.d_text_spce,'chemical shift [ppm]','normal');
    
    if strcmp(param.mode,'peak')
        str='rel. h: ';
    elseif strcmp(param.mode,'multiplet');
        str='rel. &#35;H: ';
    end
    if ts.is2c
        pd.ssHeightLegend = [...
            str,...
            '0&#60;',...
            '<tspan dy="-2.2" fill="',colhex.green.liteBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="+2.4" dx="-6" fill="',colhex.purple.liteBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="-0.4">&#60;0.2&#60;</tspan>',...
            '<tspan dy="-2.2" fill="',colhex.green.middBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="+2.4" dx="-6" fill="',colhex.purple.middBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="-0.4">&#60;0.8&#60;</tspan>',...
            '<tspan dy="-2.2" fill="',colhex.green.darkBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="+2.4" dx="-6" fill="',colhex.purple.darkBrewer,  '">&#9632;</tspan>',...
            '<tspan dy="-0.4">&#60;1</tspan>'];
    else
        pd.ssHeightLegend = [...
            str,...
            '0&#60;',...
            '<tspan dy="-0.5" dx="-0.3" fill="',colhex.green.liteBrewer,'">&#9632;</tspan>',...
            '<tspan dy="+0.5" dx="-0.3">&#60;0.2&#60;</tspan>',...
            '<tspan dy="-0.5" dx="-0.3" fill="',colhex.green.middBrewer,'">&#9632;</tspan>',...
            '<tspan dy="+0.5" dx="-0.3">&#60;0.8&#60;</tspan>',...
            '<tspan dy="-0.5" dx="-0.3" fill="',colhex.green.darkBrewer,'">&#9632;</tspan>',...
            '<tspan dy="+0.5" dx="-0.3">&#60;1</tspan>'];
    end
    svgo_text_c(3+sh2x(min(pd.xRange)),pd.d1+1.9*pd.d_text_hght+2*pd.d_text_spce,pd.ssHeightLegend,'legend','end');
    fp('</g>');
    fp('</g>');
    % ##### METABOLITE DETAILS TITLE #####
    pd.xColCas = 7*pd.d_text_wdth;
    ggo(pd.o_metabo_title);
    if ismember(ps.param.variant,{'2c','pm2c'});
        A=max(sum(cellfun(@length,spectrum.met_cas),2));
        pd.xColName = pd.xColCas+(A+5)*pd.d_text_wdth;
    else
        pd.xColName = pd.xColCas+12*pd.d_text_wdth;
    end
    svgo_text_c(0,pd.d_text_hght,'score','normal','s');
    
    pd.xColName=pd.xColCas-8;
    if nVar==2
        svgo_text_c(pd.xColName,pd.d_text_hght,[...
            'name (',...
            '<tspan dy="-0.8" fill="',colhex.green.darkBrewer,'">&#9632;</tspan>',...
            '<tspan dy=" 0.8"> &#38; </tspan>',...
            '<tspan dy="-0.8" fill="',colhex.purple.darkBrewer,  '">&#9632;</tspan>',...
            '<tspan dy=" 0.8">)</tspan>'],...
            'normal','s');
    else
        svgo_text_c(pd.xColName,pd.d_text_hght,'name','normal','s');
    end
    fp('</g>');
    % ##### METABOLITE DETAILS #####
    ggo(pd.o_metabo_boxxx);
    if ts.howto
        svgo_rect([29,130],pd.y_spectr_lines(1)-pd.d_text_line/2+pd.d_text_hght*7/9*[-1,1],'#ffffff','#e31a1c');
        svgo_text_c(132,pd.y_spectr_lines(1)-2,'F','howto','s');
    end
    Z=nshow/4;
    for i=1:(Z-1)
        svgo_line([0,pd.s_metabo_boxxx(1)+20],mean(pd.y_spectr_lines([4*i,4*i+1]))-pd.d_text_line/2,pd.c_grid,pd.dl_grid);
    end
    svgo_line([0,pd.s_metabo_boxxx(1)+20],0,pd.c_grid,pd.dl_grid);
    svgo_line([0,pd.s_metabo_boxxx(1)+20],pd.s_spectr_boxxx(2),pd.c_grid,pd.dl_grid);
    
    casExponent=0;
    for i = 1:nshow
        str=spectrum.met_cas{i};
        A=regexp(str,'-','split');
        A=floor(log10(str2double(A{1})));
        casExponent=max([casExponent,A]);
    end
    lgd={};
    for i = 1:nshow
        if max(spectrum.score)<100
            str = num2str(spectrum.score(i),'%.1f');
            if max(spectrum.score)>=10
                if spectrum.score(i)<10
                    str=add0(str);
                end
            end
        else
            str = num2str(round(spectrum.score(i)));
            if spectrum.score(i)<100
                str = add0(str);
            end
        end
        if isfield(ps,'op')
            if j==1
                fprintf('.. %s - sco %.1f - % op %.1\n',ps.tag{jPseudo},spectrum.score(i),ps.op(jPseudo));
            end
            if spectrum.score(i)>ps.op(jPseudo)
                svgo_text_c(0,pd.y_spectr_lines(i)-pd.d_text_midd,str,'permsig','s');
            else
                svgo_text_c(0,pd.y_spectr_lines(i)-pd.d_text_midd,str,'normal','s');
            end
        else
            svgo_text_c(0,pd.y_spectr_lines(i)-pd.d_text_midd,str,'normal','s');
        end
        if ts.is2c
            if length(spectrum.met_name{i,1})<14
                str1 = spectrum.met_name{i,1};
            else
                str1 = spectrum.met_cas{i,1};
            end
            if length(spectrum.met_name{i,2})<14
                str2 = spectrum.met_name{i,2};
            else
                str2 = spectrum.met_cas{i,2};
            end
            str = [str1,' &#38; ',str2];
        else
            str=spectrum.met_cas{i};
            A=regexp(str,'-','split');
            A=floor(log10(str2double(A{1})));
            while A<casExponent
                str=add0(str);
                A=A+1;
            end
        end
        if ts.is2c
            if length(spectrum.met_name{i,1})<pd.n_maxchar
                str1 = spectrum.met_name{i,1};
            else
                str1 = spectrum.met_cas{i,1};
                lgd = [lgd;{spectrum.met_cas{i,1},spectrum.met_name{i,1}}];
            end
            if length(spectrum.met_name{i,2})<pd.n_maxchar
                str2 = spectrum.met_name{i,2};
            else
                str2 = spectrum.met_cas{i,2};
                lgd = [lgd;{spectrum.met_cas{i,2},spectrum.met_name{i,2}}];
            end
            str = [str1,' &#38; ',str2];
        else
            len = length(spectrum.met_name{i});
            if ~isempty(strfind(spectrum.met_name{i},';'))
                len=len-5;
            end
            if len>25
                %                 fii=fopen('toolong','a');
                %                 fprintf(fii,'%s\n',spectrum.met_name{i});
                %                 fclose(fii);
                str = spectrum.met_cas{i};
                lgd=[lgd;{spectrum.met_cas{i},spectrum.met_name{i}}];
            else
                str = spectrum.met_name{i};
            end
        end
        if spectrum.ctrl(i)==1
            svgo_text_c(pd.xColName,pd.y_spectr_lines(i)-pd.d_text_midd,[str,'<tspan dy="-2.5" font-size="7">#</tspan>'],'control','s');
        else
            svgo_text_c(pd.xColName,pd.y_spectr_lines(i)-pd.d_text_midd,str,'normal','s');
        end
    end
    
    fp('</g>');
    
    if ~isempty(lgd)
        oldSize=ns(pd.size(2));
        [~,a]=unique(lgd(:,1));
        lgd = lgd(a,:);
        nl = ceil(length(a)/3);
        ggo([0,pd.size(2)]);
        pd.size(2)=pd.size(2)+2+nl*12+4;
        newSize=ns(pd.size(2));
        svgo_rect([0,pd.size(1)],[0,pd.size(2)],'#ffffff');
        svgo_line([0,pd.size(1)],1,pd.c_grid,pd.dl_grid);
        for i = 1:length(a)
            lv = floor((i-1)/3);
            yv = pd.size(1)*(1/6+((i-3*lv)-1)*1/3);
            svgo_text_c(yv,2+(1+lv)*pd.d_text_line-1,[lgd{i,1},': ',lgd{i,2}],'normal');
        end
        fp('</g>');
        fp('</svg>');
        fclose('all');
        file_id = fopen(fullfile(param.dir_source,[param.tag,'.svg']));
        Q=textscan(file_id,'%s','delimiter','?');
        Q=Q{1};
        Q{1}=strrep(Q{1},oldSize,newSize);
        Q{2}=strrep(Q{2},oldSize,newSize);
        fclose(file_id);
        file_id = fopen(fullfile(param.dir_source,[param.tag,'.svg']),'w');
        fprintf(file_id,'%s\n',Q{:});
        fclose(file_id);
    else
        fp('</svg>');
        fclose(file_id);
    end
end
