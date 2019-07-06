clear all;
% SET FUNCTION PATH
funcdir=getenv('DR_METABOMATCHING');
if ~isempty(funcdir)
    addpath(fullfile(funcdir,'func'));
else
    addpath('func');
end
ini__;
dirs = dir;
dirs_source = {dirs(:).name};
dirs_source = dirs_source(strncmp(dirs_source,'ps.',3));
for i = 1:length(dirs_source)
    clear ps;
    dir_source = dirs_source{i};
    fprintf('+---------------------------------------\n')
    fprintf('|   %s\n+--\n',dir_source);
    ps = fun_load_parameters(dir_source);
    fprintf('--- loading data ------------------');
    ps = fun_build_spectrum_database(ps);
    ps = fun_import_correlation(ps);
    ps = fun_import_pseudospectra(ps);
    ps = fun_processps(ps);
    fprintf(' done\n--- building transition matrix ----');
    ps = fun_build_mim(ps);    
    fprintf(' done\n--- metabomatching ----------------');
    ps = fun_metabomatching_core(ps);
    if ps.param.n_permutation>0
        fprintf(' done\n--- permutation -------------------');
        ps = fun_pshuffle(ps);
    end
    fprintf(' done\n--- writing scores ----------------');
    fun_write_scores(ps);
    save(fullfile(dir_source,'ps.mat'),'ps');
    fprintf(' done\n--- writing svg -------------------');
    vis_metabomatching(dir_source);
    fprintf(' done\n----------------------------------- ---+');
    fprintf('\n\n');
end
