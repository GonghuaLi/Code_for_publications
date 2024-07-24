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

id_atp = findRxnIDs(Recon3,'DM_atp_c_');
id_biomass = findRxnIDs(Recon3,'biomass_reaction');
option = [id_atp,0.9;id_biomass,0.9];

expr_infor = file2cell(pars.expressionFile,'\t');
%modelnames = expr_infor(1,2:end);

expressionData.Locus = expr_infor(2:end,1);
expr = cell2float(expr_infor(2:end,2:end));

modelnames = sampleInfo(:,1);
Recon3.lb(Recon3.lb<-1e-4) = -1000;
Recon3.ub(Recon3.ub>1e-4) = 1000;
Recon3.geneSymbol= upper(Recon3.geneSymbol);
for k = 1:N
    k
    outmodel = Recon3;
    %outmodel.lb = lb(:,k);
    %outmodel.ub = ub(:,k);
    exchangeIdx = sum(abs(outmodel.S))==1 & sum(outmodel.S)==-1;
    outmodel.lb(exchangeIdx)= 0;
    for i = 1:length(bloodrxns)
        outmodel.lb(findRxnIDs(outmodel,bloodrxns{i})) = exlb(i);
        outmodel.ub(findRxnIDs(outmodel,bloodrxns{i})) = exub(i);
    end
    
    expressionData.Data = expr(:,k);
    %outmodel = createTissueSpecificModel_recon2(outmodel,expressionData,1,0,[],'GIMME',option);
     [tissueModel,Rxns] = createTissueSpecificModel_recon2(outmodel,expressionData,1,0,[],'GIMME',[id_atp,0.9]);
    [tissueModel2,Rxns2] = createTissueSpecificModel_recon2(outmodel,expressionData,1,0,[],'GIMME',[id_biomass,0.9]);

    Allrxns = unique([tissueModel.rxns;tissueModel2.rxns]);
    rmRxnsID = ~ismember(Recon3.rxns,Allrxns);

    outmodel = removeRxns(outmodel,outmodel.rxns(rmRxnsID));
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



