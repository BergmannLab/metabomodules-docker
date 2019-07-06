function ps=fun_import_pseudospectra(ps)
% FUN_IMPORT_PSEUDOSPECTRA Read pseudospectrum files
global sFLD
if exist(ps.param.dir_source,'dir')
    fn=dir(fullfile(ps.param.dir_source,'*.pseudospectrum.tsv'));
    if ~isempty(fn)
        for f=sFLD
            if isfield(ps,f{1});
                ps = rmfield(ps,f{1});
            end
        end
        for i = 1:length(fn)
            dsfni = fullfile(ps.param.dir_source,fn(i).name);            
            ts = fun_readtable(dsfni);
% a "/" indicates that the pseudospectrum.tsv file contains 
% multiple pseudospectra; the header is then of the form "variable/tag"            
            if ~ismember('/',[ts.lb{:}])               
                ps.tag{i,1} = strrep(fn(i).name,'.pseudospectrum.tsv','');
                for i_field=1:length(ts.lb)
                    ps.(ts.lb{i_field})(:,i)=ts.pr{i_field};
                end
            else
                ps.tag={};
                ps.param.multi=fn(i).name;
                for i_field = 1:length(ts.lb)
                    if strcmpi(ts.lb{i_field},'shift')
                        ps.shift=ts.pr{i_field};
                    else
                        xx=regexp(ts.lb{i_field},'/','split');
                        if length(xx)==2
                            ix = find(strcmpi(ps.tag,xx{2}));
                            if isempty(ix)
                                ix = length(ps.tag)+1;
                                ps.tag{ix,1} = xx{2};
                            end
                            ps.(xx{1})(:,ix)=ts.pr{i_field};
                        end
                    end
                end
            end
        end
        ps.shift = ps.shift(:,1);
        %
        if ismember('beta',fieldnames(ps))
            ps.param.pstype='xs';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='p';
            end
%             if ~isfield(ps.param,'significant')
%                 ps.param.significant=-log10(5E-8);
%             else
%                 if ps.param.significant<1
%                     ps.param.significant=-log10(ps.param.significant);
%                 end
%             end
        elseif ismember('isa',fieldnames(ps)); 
            ps.param.pstype='isa';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='z';
            end
        elseif ismember('cr',fieldnames(ps));
            if ~ismember('samplesize',fieldnames(ps.param))
                % warning('metabomatching:needSS','For correlation-type pseudospectra, the sample size needs to be passed as a parameter [''samplesize''] to metabomatching for proper scaling. Without the sample size, the Fisher transformed correlations are simply standardized.');
            end
            ps.param.pstype='correlation';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='z';
            end
        elseif ismember('pc',fieldnames(ps))
            ps.pca=ps.pc;
            ps=rmfield(ps,'pc');
            ps.param.pstype='pca';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='z';
            end
        elseif ismember('pca',fieldnames(ps));
            ps.param.pstype='pca';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='z';
            end        
        elseif ismember('z',fieldnames(ps));
            ps.param.pstype='z';
            if ~isfield(ps.param,'plot_type')
                ps.param.plot_type='z';
            end        
        else          
            error('metabomatching:noTYPE','Couldn''t identify the pseudospectrum type.');
        end        
    else
        error('metabomatching:noPS','No pseudospectrum files found (*.pseudospectrum.tsv)');
    end
else
    error('metabomatching:noDS','Source directory %s does not exist',ps.param.dir_source);
end
