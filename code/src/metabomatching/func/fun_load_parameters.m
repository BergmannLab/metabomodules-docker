function ps = fun_load_parameters(dir_source,fn)
% FUNCTION_LOAD_PARAMETERS Read parameter file
ps.param.dir_source = dir_source;
if nargin<2
    fn = fullfile(ps.param.dir_source,'parameters.in.tsv');
end
global sNF
number_fields=[sNF,'p_significant','p_suggestive'];
defaults = { ...
    'variant','1c';...
    'mode','peak';...
    'scoring','chisq';...
    'reference','hmdb';...
    'decorr_lambda',1;...
    'n_show',8;...
    'n_permutation',0};
%
special_defaults.dsh_peak = 0.025;
special_defaults.dsh_multiplet = 0.010;
special_defaults.suggestive = 5;
% /
if exist(fn,'file');
    pr = fun_read(fn,'%s%s');
    for j = 1:length(pr{1})
        field = pr{1}{j};
        value = pr{2}{j};
        if ismember(field,number_fields)
            ps.param.(field)=str2double(value);
        else
            field=strrep(field,'.','_');
            ps.param.(field)=value;
        end
    end
else
    fprintf('|   no parameter file, using defaults\n---\n');
end
%
%
% ##### PARAMETERS THAT CHANGED NAME #####
if isfield(ps.param,'p_suggestive')
  ps.param.suggestive = ps.param.p_suggestive;
  ps.param = rmfield(ps.param,'p_suggestive');
end
if isfield(ps.param,'p_significant')
  ps.param.significant = ps.param.p_significant;
  % warning('The parameter ''p_significant'' has been replaced by ''significant''.');
  ps.param = rmfield(ps.param,'p_significant');
end
if isfield(ps.param,'significant');
    if ps.param.significant<1
        ps.param.significant=-log10(ps.param.significant);
    end
end
%
%
% ##### ASSIGN DEFAULTS #####
for i = 1:size(defaults,1)
    if ~isfield(ps.param,defaults{i,1})
        ps.param.(defaults{i,1}) = defaults{i,2};
    end
end
if ~isfield(ps.param,'dsh')
    if strcmp(ps.param.mode,'peak');
        ps.param.dsh=special_defaults.dsh_peak;
    else
        ps.param.dsh=special_defaults.dsh_multiplet;
    end
end
if ismember(ps.param.variant,{'pm','pm1c','pm2c'}) && ...
        ~isfield(ps.param,'suggestive')
    ps.param.suggestive=special_defaults.suggestive;
end
ps.param.dir_source = dir_source;
