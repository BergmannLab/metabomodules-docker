function ps = fun_metabomatching_core(ps)
% FUNCTION_METABOMATCHING_CORE  Compute metabomatching scores
nm = size(ps.sid,1);
ns = size(ps.z,2);
ps.score = NaN(nm,ns);
L = ps.param.decorr_lambda;
for jm = 1:nm
    if L<1
        se_met=[];
        for ii = 1:size(ps.cluster{jm},1)
            se_clust = find(...
                ps.shift<=ps.cluster{jm}(ii,2) & ...
                ps.shift>=ps.cluster{jm}(ii,1));
            subC{ii} = ps.correlation(se_clust,se_clust);
            se_met = [se_met;se_clust];
        end
        nb_inmet = length(se_met);
        D = blkdiag(subC{:});
        clear subC
        C = (1-L)*D+L*eye(nb_inmet);
        z = ps.z(se_met,:);
        switch ps.param.scoring
            case 'chisq'
                chisq = sum(z.*(C\z),1);
                ps.score(jm,:) = -log10(gamcdf_tail(chisq,nb_inmet/2,2));
                se = ps.score(jm,:)==Inf;
                if ~isempty(se)
                    ps.score(jm,se)=gamcdf_bound(chisq(se),nb_inmet);
                end
            case 'z'
                ps.score(jm,:) = -log10(2*normcdf(-abs(sum(z,1)),0,sqrt(abs(sum(C(:))))));
        end
    else
        se_met = find(ps.mim(:,jm));
        nb_inmet = length(se_met);
        z = ps.z(se_met,:);
        switch ps.param.scoring
            case 'chisq'
                chisq = sum(z.*z,1);
                ps.score(jm,:) = -log10(gamcdf_tail(chisq,nb_inmet/2,2));
                se = find(ps.score(jm,:)==Inf);
                if ~isempty(se)
                    ps.score(jm,se)=gamcdf_bound(chisq(se),nb_inmet);
                end
            case 'z'
                ps.score(jm,:) = -log10(2*normcdf(-abs(sum(z,1)),0,sqrt(nb_inmet)));
        end
    end
end
