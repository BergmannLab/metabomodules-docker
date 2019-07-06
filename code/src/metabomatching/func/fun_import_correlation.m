function ps = fun_import_correlation(ps)
% FUNCTION_IMPORT_CORRELATION  Import feature-feature correlation matrix

if ps.param.decorr_lambda<1
    fn = fullfile(ps.param.dir_source,'correlation.csv');
    if exist(fn,'file')
        ps.correlation=csvread(fullfile(ps.param.dir_source,'correlation.csv'));        
    else
        ps.param.decorr_lambda=1;
        ps.param.correlation_missing=1;
    end
end
