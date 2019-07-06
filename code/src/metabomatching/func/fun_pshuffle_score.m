function ps = fun_pshuffle_score(ps)
% FUNCTION_PPERMUTE
nm = length(ps.sid);
nf = length(ps.shift);
ns = length(ps.tag);
np = ps.param.n_permutation;

for js = 1:ns
  
  permps = ps;
  permps.z = squeeze(ps.zperm(:,:,js));
  
  permps = fun_metabomatching_core(permps);
  opi = permps.score;

  % Set NaN or Inf cases to 0; this will disregard them for the permscore.
  opi(~isfinite(opi))=0;
  opi = sort(opi,1,'descend');
  ps.permscore(:,js)=max(opi,[],1);
  tm = sum(repmat(ps.score(:,js),1,np)<repmat(ps.permscore(:,js)',nm,1),2);
  ps.scoreadj(:,js)=abs(log10((1+tm)/(1+np)));
end