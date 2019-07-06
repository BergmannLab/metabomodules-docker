function ps = fun_build_mim(ps)
% FUN_BUILD_MIM  Build metabolite match set matrix
fn = fullfile(ps.param.dir_source,'metdb.mat');
load(fn);

ps.sid = metdb.sid;
ps.pos = metdb.pos;
nf = length(ps.shift);
typ = ps.param.mode;
if ismember(ps.param.variant,{'2c','pm2c'});
    typ=[typ,'-2c'];
end
switch typ
    case 'peak'
        nm = length(ps.sid);
        ps.mim = zeros(nf,nm);
        for jm = 1:nm
            [tt.pos{jm},reo] = sort(ps.pos{jm});
            jCluster = 0;
            prvShift = -100;
            for jj = 1:length(tt.pos{jm})
                if tt.pos{jm}(jj)-prvShift <= 2*ps.param.dsh
                    ps.cluster{jm,1}(jCluster,2) = tt.pos{jm}(jj)+ps.param.dsh;
                else
                    jCluster = jCluster+1;
                    ps.cluster{jm,1}(jCluster,:) = tt.pos{jm}(jj)+ps.param.dsh*[-1,1];
                end
                prvShift = tt.pos{jm}(jj);
                sl = abs(ps.shift-tt.pos{jm}(jj))<=ps.param.dsh;
                ps.mim(sl,jm) = 1;
            end
        end
    case 'multiplet'
        nm = length(ps.sid);
        ps.mim = zeros(nf,nm);
        for jm = 1:nm
            ps.cluster{jm} = ps.pos{jm}+ps.param.dsh*repmat([-1,1],size(ps.pos{jm},1),1);
            [~,reo] = sort(ps.cluster{jm}(:,1));
            ps.cluster{jm}=ps.cluster{jm}(reo,:);
            for jj = 1:size(ps.pos{jm},1)
                if jj>1
                    A=ps.cluster{jm}(jj-1,2);
                    B=ps.cluster{jm}(jj  ,1);
                    if A>=B
                        ps.cluster{jm}(jj-1,2)=mean([A,B])-1e-5;
                        ps.cluster{jm}(jj  ,1)=mean([A,B]);
                    end
                end
                sl1 = ps.shift>=(ps.pos{jm}(jj,1)-ps.param.dsh);
                sl2 = ps.shift<=(ps.pos{jm}(jj,2)+ps.param.dsh);
                ps.mim(sl1&sl2,jm) = 1;
            end
        end
    case 'peak-2c'
        tt.sid = ps.sid;
        %tt.sid = ps.sid;
        %tt.idDat = ps.sidDat;
        tt.pos = ps.pos;
        nt = length(tt.sid);
        tt.mim = zeros(nf,nt);
        for jm = 1:nt
            for jj = 1:length(tt.pos{jm})
                sl = abs(ps.shift-tt.pos{jm}(jj))<=ps.param.dsh;
                tt.mim(sl,jm) = 1;
            end
        end
        nm = nt*(nt+1)/2;
        ps.is = zeros(nm,2);
        ps.pos = cell(nm,2);
        ps.mim = zeros(nf,nm);
        k = 0;
        for jm1 = 1:nt
            tic;
            st = jm1:nt;
            n0 = length(st);
            tt.mim1 = repmat(tt.mim(:,jm1),1,n0);
            tt.mim2 = tt.mim(:,st);
            ps.mim(:,(k+1):(k+n0))=tt.mim1|tt.mim2;
            for F={'sid','pos'}%,'sid','idDat'};
                f=F{1};
                ps.(f)((k+1):(k+n0),1) = repmat(tt.(f)(jm1),n0,1);
                ps.(f)((k+1):(k+n0),2) = tt.(f)(st);
            end
            k=k+length(st);
        end
        for jm = 1:nm
            [tt,reo] = sort([ps.pos{jm,1};ps.pos{jm,2}]);
            jCluster = 0;
            prvShift = -100;
            for jj = 1:length(tt)
                if tt(jj)-prvShift <= 2*ps.param.dsh
                    ps.cluster{jm,1}(jCluster,2) = tt(jj)+ps.param.dsh;
                else
                    jCluster = jCluster+1;
                    ps.cluster{jm,1}(jCluster,:) = tt(jj)+ps.param.dsh*[-1,1];
                end
                prvShift = tt(jj);
            end
        end
%         [ism,loc]=ismember(ps.sid(:,1),metdb.sid);
%         ps.cas=cell(nm,2);
%         ps.sid=zeros(nm,2);
%         ps.cas(ism,1)=metdb.cas(loc(ism));
%         ps.sid(ism,1)=metdb.sid(loc(ism));
%         [ism,loc]=ismember(ps.sid(:,2),metdb.sid);
%         ps.cas(ism,2)=metdb.cas(loc(ism));
%         ps.sid(ism,2)=metdb.sid(loc(ism));
    case 'multiplet-2c'
        tt.sid = ps.sid;
        tt.pos = ps.pos;
        nt = length(tt.sid);
        tt.mim = zeros(nf,nt);
        for jm = 1:nt
            for jj = 1:size(tt.pos{jm},1)
                sl1 = ps.shift>=(tt.pos{jm}(jj,1)-ps.param.dsh);
                sl2 = ps.shift<=(tt.pos{jm}(jj,2)+ps.param.dsh);
                tt.mim(sl1&sl2,jm) = 1;
            end
        end
        nm = nt*(nt+1)/2;
        ps.sid = zeros(nm,2);
        ps.mim = zeros(nf,nm);
        k = 0;
        for jm1 = 1:nt
            tic;
            st = jm1:nt;
            n0 = length(st);
            tt.mim1 = repmat(tt.mim(:,jm1),1,n0);
            tt.mim2 = tt.mim(:,st);
            ps.mim(:,(k+1):(k+n0))=tt.mim1|tt.mim2;
            ps.sid((k+1):(k+n0),1)=repmat(tt.sid(jm1),n0,1);
            ps.sid((k+1):(k+n0),2)=tt.sid(st);
            
            for ijm = 1:n0
                jm2 = st(ijm);
                jCluster=0;
                clust = 0;
                for jj = 1:size(ps.mim,1)
                    if ps.mim(jj,k+ijm)==1 && clust==0
                        jCluster=jCluster+1;
                        ps.cluster{k+ijm,1}(jCluster,1)=ps.shift(jj);
                        clust=1;
                    end
                    if ps.mim(jj,k+ijm)==0 && clust==1
                        ps.cluster{k+ijm,1}(jCluster,2)=ps.shift(jj-1);
                        clust=0;
                    end
                end
            end
            k=k+length(st);
        end

end

[nm,nk]=size(ps.sid);
ps.cas=cell(nm,nk);
for ik = 1:nk
        [ism,loc]=ismember(ps.sid(:,ik),metdb.sid);
        ps.cas(ism,ik)=metdb.cas(loc(ism));
end


