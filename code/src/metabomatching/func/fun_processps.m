function ps = fun_processps(ps)
global sFLD sNF
% defining cut points when computing permutation scores requires sorted shifts.
[~,r] = sort(ps.shift);
sfld = [sFLD,'shift'];
for jd = 1:length(sfld)
    fld = sfld{jd};
    if ismember(fld,fieldnames(ps))
        ps.(fld)=ps.(fld)(r,:);
    end
end

switch ps.param.pstype
    case 'xs'
        % remove dummy columns
        se = find(all(ps.beta==0)|all(ps.p==0)|all(ps.se)==0);
        ps.tag(se)=[];
        ps.beta(:,se)=[];
        ps.se(:,se)=[];
        ps.p(:,se)=[];
        % metabomatching only uses beta and x in their ratio.
        ps.z=ps.beta./ps.se;
    case 'isa'
        se = find(all(ps.isa==0));
        ps.tag(se)=[];
        ps.isa(:,se)=[];
        ps.z=zscore(ps.isa);
    case 'pca'
        se = find(all(ps.pca==0));
        ps.tag(se)=[];
        ps.pca(:,se)=[];
        ps.z=zscore(ps.pca);
    case 'correlation'
        se = find(all(ps.cr==0));
        ps.tag(se)=[];
        ps.cr(:,se)=[];
        if isfield(ps.param,'samplesize')
          ps.z = atanh(ps.cr)*sqrt(ps.param.samplesize-3);
          if isfield(ps.param,'crscale')
            ps.z = ps.z/ps.param.crscale;
          else
            ps.z = zscore(ps.z);
          end
        else
          ps.z = zscore(atanh(ps.cr));
        end
end

% handling pm cases
if ismember(ps.param.variant,{'pm','pm1c','pm2c'})
  ps.param.pmkeeporig=true;
  F = fieldnames(ps);
  nr = length(ps.tag);
  nosig = ~isfield(ps.param,'significant');
  nosug = ~isfield(ps.param,'suggestive');
  for jr = 1:nr
        
    if nosig
      sig = max(4,prctile(abs(ps.z(:,jr)),95));
    else
      sig = ps.param.significant;
    end
    
    if nosug
      sug = max(3,prctile(abs(ps.z(:,jr)),85));
    else
      sug = ps.param.suggestive;
    end
    
    % 
	  if strcmp(ps.param.pstype,'xs')
      se_sug_neg = ps.p(:,jr)<10^(-sug) & ps.beta(:,jr)<0;
      se_sug_pos = ps.p(:,jr)<10^(-sug) & ps.beta(:,jr)>0;
      sig_pos = any(ps.p(:,jr)<10^(-sig) & ps.beta(:,jr)>0);
      sig_neg = any(ps.p(:,jr)<10^(-sig) & ps.beta(:,jr)<0);
	  else
  		se_sug_neg = ps.z(:,jr)<-sug;
      se_sug_pos = ps.z(:,jr)>+sug;
      sig_pos = any(ps.z(:,jr)>+sig);
      sig_neg = any(ps.z(:,jr)<-sig);
	  end
    if sig_pos && sig_neg
      
      FF = fieldnames(ps);
      for jf = 1:length(FF)
        ff = FF{jf};
        if ~ismember(ff,{'shift','param','tag'})
          ps.(ff)=[ps.(ff),ps.(ff)(:,jr)];          
          ps.(ff)=[ps.(ff),ps.(ff)(:,jr)];
        end
      end
      ps.tag = [ps.tag;[ps.tag{jr},'.neg']];
      if ps.param.pmkeeporig
        ps.tag = [ps.tag;[ps.tag{jr},'.pos']];
      else
        ps.tag{jr} = [ps.tag{jr},'.pos'];
      end
      if ps.param.pmkeeporig
        ps.z(se_sug_neg,end  )=-1E-10;
        ps.z(se_sug_pos,end-1)=+1E-10;
      else
        ps.z(se_sug_neg,jr )=-1E-10;
        ps.z(se_sug_pos,end)=+1E-10;
      end
    
      % this may not be necessary, even for xs-type
		  if strcmp(ps.param.pstype,'xs')
  	    ps.beta(se_sug_neg,jr) =-1E-10;
        ps.beta(se_sug_pos,end)=+1E-10;
	    end
    end
  end
end
