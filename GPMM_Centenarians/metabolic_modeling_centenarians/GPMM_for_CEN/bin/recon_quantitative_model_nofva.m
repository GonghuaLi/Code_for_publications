function [modelout,parout] = recon_quantitative_model_nofva(annoted_Recon3_model,enzyme_expr,varargin )
%% reconstruction of quantitiative model by protein aboundance

modelout = annoted_Recon3_model;
fluxByBlood =  [modelout.lb,modelout.ub];
%% check if essentail genes are expressed
IsEssentail = annoted_Recon3_model.IsEssentail;
idx = enzyme_expr<1 & IsEssentail>0;
% Essentail analysis
if sum(idx)>=1
    enzyme_expr(idx) = 1.001;
end

%%  constraint by enzym activity
rxns_expr = get_rxns_expr_byrule(modelout,enzyme_expr);
[IA,IB] = ismember(modelout.rxns,annoted_Recon3_model.rxns);
modelout_Kcat = annoted_Recon3_model.Kcat(IB);
rxns_active = ones(length(modelout.rxns),1)*-1000;
for i = 1:length(rxns_active)
    if modelout.IsExprPQM(i) == 0 
        rxns_active(i) = 0;
    elseif modelout.IsExprPQM(i) > 0 & (modelout_Kcat(i) <= 0 | rxns_expr(i) <= 0)
        continue;
    elseif modelout.IsExprPQM(i) > 0 & modelout_Kcat(i) * rxns_expr(i)  < 1e-9
        continue;
    else
        rxns_active(i) = modelout_Kcat(i) * rxns_expr(i);
    end
end
idx = fluxByBlood(:,2)<1000 & fluxByBlood(:,2)>0 & rxns_active >= 1e-9;
if ~isempty(varargin)
    par = varargin{1};
else
    par = median(fluxByBlood(idx,2)./rxns_active(idx));
end
parout = par;
ub1 = par*rxns_active;
ub1(ub1<0) = 1000;
lb1 = -ub1;
modelout.ub = min([ub1,fluxByBlood(:,2)],[],2);
modelout.lb = max([lb1,fluxByBlood(:,1)],[],2);

% set Oxidative phosphorylation active
%Oxidative = {'ATPS4m','CYOOm2','CYOR_u10m','NADH2_u10m',...
%    'r0205','CYOOm3','FADH2ETC','GLYC3PFADm'};
Oxidative = {'ATPS4mi','CYOOm2i','CYOR_u10mi','NADH2_u10mi',...
    'r0205','CYOOm3i','FADH2ETC','GLYC3PFADm'};

id_oxidative = findRxnIDs(modelout,Oxidative);
modelout.lb(id_oxidative) = 0;
modelout.ub(id_oxidative) = 1000;

%gg = changeObjective(modelout,'biomass_reaction');
%optimizeCbModel(gg)



