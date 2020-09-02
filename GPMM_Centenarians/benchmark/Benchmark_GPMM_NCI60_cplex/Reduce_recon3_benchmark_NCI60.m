%Reduce_recon3_benchmark_NCI60
pars = parse_parsfile('./pars.txt');
load(pars.curatedRecon3);
%%
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
%%
changeCobraSolver('cplex');
outmodel = Recon3;
idATP = findRxnIDs(outmodel,'DM_atp_c_');
idbiomass = findRxnIDs(outmodel,'biomass_reaction');
%% exchange data
vv= file2cell('./data/NCI60_flux_benchmark_Recon3(mmol L min)59.txt','\t');
bloodrxns = vv(2:end,1);
fluxes = cell2float(vv(2:end,2:end));
%%
vid = findRxnIDs(Recon3,bloodrxns) >0;
bloodrxns = bloodrxns(vid);
fluxes = fluxes(vid,:);
%% uptake 
e1 = mean(fluxes,2);
iduptake = e1 <0;
uptakeids = bloodrxns(iduptake);
exlb = e1(iduptake);
exub = repmat(1000,length(exlb),1);
addtionUptake = {'EX_his_L(e)','EX_pi(e)','EX_na1(e)','EX_fe2(e)','EX_k(e)',...
    'EX_h2o(e)','EX_Tyr_ggn(e)','EX_o2(e)'};
x1 = [-0.05;repmat(-1000,length(addtionUptake)-1,1)];
exlb= [exlb;x1];
exub= [exub;repmat(1000,length(addtionUptake),1)];
uptakeids = [uptakeids; addtionUptake'];   
%%
exchangeIdx = sum(abs(outmodel.S))==1 & sum(outmodel.S)==-1;
outmodel.lb(exchangeIdx)= 0;
for i = 1:length(uptakeids)
   % bloodrxns{i}
    outmodel.lb(findRxnIDs(outmodel,uptakeids{i})) = exlb(i);
    outmodel.ub(findRxnIDs(outmodel,uptakeids{i})) = exub(i);
end
%%
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
flux = FastMM_FVA(mm);
rmrxns = mm.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
gg = removeRxns(mm,rmrxns);

%% write new recon 3
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










