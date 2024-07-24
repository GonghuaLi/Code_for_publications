pars = parse_parsfile('pars.txt');

changeCobraSolver(pars.cobrasolver)
load(pars.curatedRecon3)
load('./data/lb_ub_for_PQMM');
vv = file2cell(pars.exchangeFile,'\t');
%%
id = regexp(vv(:,7),'^EX_');
id = ~isemptycell(id);
bloodrxns = vv(id,7);
exlb = cell2float(vv(id,11));
exub = cell2float(vv(id,12));
N = size(lb,2);

if ~isdir('./Mat_PQM')
    mkdir('./Mat_PQM')
end

modelnames = sampleInfo(:,1);
for k = 1:N
    k
    outmodel = Recon3;
    outmodel.lb = lb(:,k);
    outmodel.ub = ub(:,k);
    exchangeIdx = sum(abs(outmodel.S))==1 & sum(outmodel.S)==-1;
    outmodel.lb(exchangeIdx)= 0;
    for i = 1:length(bloodrxns)
        outmodel.lb(findRxnIDs(outmodel,bloodrxns{i})) = exlb(i);
        outmodel.ub(findRxnIDs(outmodel,bloodrxns{i})) = exub(i);
    end
    %if strcmp(pars.degradation,'yes')
    %    outmodel.ub(findRxnIDs(outmodel,'r1333')) = 1000;
    %else
    %    outmodel.ub(findRxnIDs(outmodel,'r1333')) = 0;
    %end
    if strcmp(pars.isFastMM,'yes')
        
        flux = FastMM_FVA(outmodel);
    else
        [minf, maxf] = fluxVariability(outmodel);
	flux = [minf,maxf];
    end
    %optimizeCbModel(outmodel)
    %flux = FastMM_FVA_gurobi5(outmodel);
    outmodel.lb = flux(:,1);
    outmodel.ub = flux(:,2);
    rmRxns = outmodel.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
    outmodel = removeRxns(outmodel,rmRxns);
    %optimizeCbModel(outmodel)
    %stat = cobra2FastKO(outmodel,['./FastKO_PQM_add/',modelnames{k}]);
    save(['./Mat_PQM/',modelnames{k},'.mat'],'outmodel');
end



