function outmodel = Reduce_recon3(parsfile)
pars = parse_parsfile(parsfile);
load(pars.curatedRecon3);
idTrans = regexp(Recon3.subSystems,'Transport');
idExchange = regexp(Recon3.subSystems,'Exchange');
idOxidation = regexp(Recon3.subSystems,'Oxidative');
idTrans = ~isemptycell(idTrans);
idExchange = ~isemptycell(idExchange);
idOxidation = ~isemptycell(idOxidation);
ExpressedID = Recon3.Kcat>0 | idTrans' | idExchange' | idOxidation';
UnExpressedID = ~ExpressedID;
ExpressedRxns = Recon3.rxns(ExpressedID);
UnExpressedRxns = Recon3.rxns(UnExpressedID);
%%  add tissue specific objective function ----important
changeCobraSolver(pars.cobrasolver);
outmodel = Recon3;
idATP = findRxnIDs(outmodel,'DM_atp_c_');
idbiomass = findRxnIDs(outmodel,'biomass_reaction');
%% read exchange data
vv = file2cell(pars.exchangeFile,'\t'); % could be change
id = regexp(vv(:,7),'^EX_');
id = ~isemptycell(id);
bloodrxns = vv(id,7);
exlb = cell2float(vv(id,11));
exub = cell2float(vv(id,12));
exchangeIdx = sum(abs(outmodel.S))==1 & sum(outmodel.S)==-1;
outmodel.lb(exchangeIdx)= 0;
for i = 1:length(bloodrxns)
   % bloodrxns{i}
    outmodel.lb(findRxnIDs(outmodel,bloodrxns{i})) = exlb(i);
    outmodel.ub(findRxnIDs(outmodel,bloodrxns{i})) = exub(i);
end


%%
%outmodel.ub(findRxnIDs(outmodel,'r1333')) = 0;
%[tissueModel,Rxns] = gimme_by_reaction(outmodel, ExpressedRxns,UnExpressedRxns,[],[idATP,0.9;idbiomass,0.9]);
[tissueModel,Rxns] = gimme_by_reaction(outmodel, ExpressedRxns,UnExpressedRxns,[],[idATP,0.9]);
[tissueModel2,Rxns2] = gimme_by_reaction(outmodel, ExpressedRxns,UnExpressedRxns,[],[idbiomass,0.9]);
%%
Allrxns = unique([tissueModel.rxns;tissueModel2.rxns]);
rmRxnsID = ~ismember(Recon3.rxns,Allrxns);

mm = removeRxns(outmodel,outmodel.rxns(rmRxnsID));
mm = changeObjective(mm,'DM_atp_c_');
optimizeCbModel(mm)
mm2 = changeObjective(mm,'biomass_reaction');
optimizeCbModel(mm2)
%%
mm3 = changeObjective(mm,'DM_atp_c_');
f = optimizeCbModel(mm3)
mm3.lb(findRxnIDs(mm3,'DM_atp_c_')) = f.f*0.9;
mm3 = changeObjective(mm3,'biomass_reaction');
optimizeCbModel(mm3)
%gg = mm;
%%
%flux = FastMM_FVA(mm);
if strcmp(pars.isFastMM,'yes')
        
    flux = FastMM_FVA(mm);
else
    [minf, maxf] = fluxVariability(mm);
    flux = [minf,maxf];
end
%flux = FastMM_FVA_gurobi5(mm);
rmrxns = mm.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
gg = removeRxns(mm,rmrxns);

%% write new recon 2
IsExprPQM = ismember(Recon3.rxns,gg.rxns);
Recon3.IsExprPQM = IsExprPQM;
gg = changeObjective(gg,'biomass_reaction');
optimizeCbModel(gg)
%koFlux = FastMM_singleGeneKO(gg);
%[xx,fmax] = singleGeneDeletion(gg);
%[tmp,koflux] = singleGeneDeletion(gg);
[grRatio, grRateKO, grRateWT]=singleGeneDeletion(gg);
maxBiomass = grRateWT;%koFlux(1,3);
Recon3.IsEssentail = grRateKO < maxBiomass * 0.6;%koFlux(2:end,3) < maxBiomass*0.6;
%Recon3.ub(findRxnIDs(Recon3,'C4STMO2r')) = 1000;%ERRO BY LGH
%ia = ismember(Recon3.subSystems,'Bile acid synthesis');%not liver delet
%Recon3.lb(ia) = 0;
%Recon3.ub(ia) = 0;

%%
%flux = FastMM_FVA(tissueModel);
%rmrxns = tissueModel.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
%Reducedmodel = removeRxns(tissueModel,rmrxns);
%% write new recon 3
%IsExprPQM = ismember(Recon3.rxns,Reducedmodel.rxns);
%Recon3.IsExprPQM = IsExprPQM;
%Recon3.IsExprPQM(findRxnIDs(outmodel,'r1333'))=1; % for burned muscle
%save(pars.curatedRecon3,'Recon3');
save(pars.curatedRecon3,'Recon3');
outmodel = Recon3;
