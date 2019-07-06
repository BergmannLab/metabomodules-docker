function fun_write_scores(ps)
% FUN_WRITE_SCORES  Write metabomatching scores to file

if ps.param.n_permutation>0
    sfld = {'score','scoreadj'};
else
    sfld = {'score'};
end

ts2C = ismember(ps.param.variant,{'pm2c','2c'});
tsMULTI = isfield(ps.param,'multi');

nm = length(ps.sid);
ns = length(ps.tag);

if tsMULTI
    for jfld = 1:length(sfld)
        fld = sfld{jfld};
        X = ps.(fld);
        fn = fullfile(ps.param.dir_source,strrep(ps.param.multi,'pseudospectrum',fld));
        fi = fopen(fn,'w');
        if ts2C
            fprintf(fi,'cas\tid\tcas\tid');
        else
            fprintf(fi,'cas\tid');
        end
        for js = 1:ns
            fprintf(fi,'\t%s/%s',fld,ps.tag{js});
        end
        fprintf(fi,'\n');
        for jm = 1:nm
            if ts2C
                fprintf(fi,'%s\t%d\t%s\t%d',ps.cas{jm,1},ps.sid(jm,1),ps.cas{jm,2},ps.sid(jm,2));
                for js = 1:ns
                    fprintf(fi,'\t%.4f',X(jm,js));
                end
                fprintf(fi,'\n');
            else
                fprintf(fi,'%s\t%d',ps.cas{jm},ps.sid(jm));
                for js = 1:ns
                    fprintf(fi,'\t%.4f',X(jm,js));
                end
                fprintf(fi,'\n');
            end
        end
        fclose(fi);
    end
else
    for jfld = 1:length(sfld)
        fld = sfld{jfld};
        X = ps.(fld);
        for js = 1:ns
            fn=fullfile(ps.param.dir_source,[ps.tag{js},'.',fld,'.tsv']);
            fi=fopen(fn,'w');
            fprintf(fi,'cas\tid\t%s\n',fld);
            for jm = 1:nm
                if ts2C
                    fprintf(fi,'%s\t%d\t%s\t%d\t%.4f\n',...
                        ps.cas{jm,1},ps.sid(jm,1),...
                        ps.cas{jm,2},ps.sid(jm,2),X(jm,js));
                else
                    fprintf(fi,'%s\t%d\t%.4f\n',ps.cas{jm},ps.sid(jm),X(jm,js));
                end
            end
            fclose(fi);
        end
    end
end
global sNF
fi=fopen(fullfile(ps.param.dir_source,'parameters.out.tsv'),'w');
sfld=fieldnames(ps.param);
for jfld = 1:length(sfld)
    fld=sfld{jfld};
    if ismember(fld,sNF)
        format = '%.4g';
    else
        format = '%s';
    end
    fprintf(fi,['%s\t',format,'\n'],fld,ps.param.(fld));
end
fclose(fi);
