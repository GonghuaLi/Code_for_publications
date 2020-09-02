% summary PQMM results
pars = parse_parsfile('pars.txt');
% get_individual_result V0.1
load(pars.curatedRecon3);
load('./data/lb_ub_for_PQMM.mat');
%changeCobraSolver('gurobi5');
changeCobraSolver(pars.cobrasolver)
[num_mets,num_rxns] = size(Recon3.S);

modelnames = sampleInfo(:,1);
N_samples = length(modelnames);

fluxRxnsMedian = zeros(num_rxns,N_samples);
fluxRxnsMean = zeros(num_rxns,N_samples);
fluxRxnsMin = zeros(num_rxns,N_samples);
fluxRxnsMax = zeros(num_rxns,N_samples);
fluxFVAMin = zeros(num_rxns,N_samples);
fluxFVAMax = zeros(num_rxns,N_samples);

consumptionFluxMetMedian = zeros(num_mets,N_samples);
consumptionFluxMetMean = zeros(num_mets,N_samples);
consumptionFluxMetMin = zeros(num_mets,N_samples);
consumptionFluxMetMax = zeros(num_mets,N_samples);
productionFluxMetMedian = zeros(num_mets,N_samples);
productionFluxMetMean = zeros(num_mets,N_samples);
productionFluxMetMin = zeros(num_mets,N_samples);
productionFluxMetMax = zeros(num_mets,N_samples);
detaFluxMetMedian = zeros(num_mets,N_samples);
detaFluxMetMean = zeros(num_mets,N_samples);
detaFluxMetMin = zeros(num_mets,N_samples);
detaFluxMetMax = zeros(num_mets,N_samples);

for i = 1:length(modelnames)
    if exist(['./Mat_PQM/MCMC',modelnames{i},'_1.mat'],'file')   
        samples = loadSamples(['./Mat_PQM/MCMC',modelnames{i}],1,2000);
    else
         continue;
    end
    load(['./Mat_PQM/',modelnames{i},'.mat']);
    rxns = outmodel.rxns;
    [IA,IB] = ismember(Recon3.rxns,rxns);
    IB(IB==0) = length(rxns)+1;
    samples = loadSamples(['./Mat_PQM/MCMC',modelnames{i}],1,2000);
    % rm all zero points%
    samples = samples(:,sum(abs(samples),1) > 1e-5);
    %
    medianFlux = median(samples,2);
    meanFlux = mean(samples,2);
    minFluxMCMC = min(samples,[],2);
    maxFluxMCMC = max(samples,[],2);
    medianFlux = [medianFlux;0];
    meanFlux = [meanFlux;0];
    minFluxMCMC = [minFluxMCMC;0];
    maxFluxMCMC = [maxFluxMCMC;0];
    %fluxMinMax = [fluxMinMax;0,0];
    
    fluxRxnsMedian(:,i) = medianFlux(IB);
    fluxRxnsMean(:,i) = meanFlux(IB);
    fluxRxnsMin(:,i) = minFluxMCMC(IB);
    fluxRxnsMax(:,i) = maxFluxMCMC(IB);
    %fluxFVAMin(:,i) = fluxMinMax(IB,1);
    %fluxFVAMax(:,i) = fluxMinMax(IB,2);
    
    %
    mets = outmodel.mets;
    [IA,IB] = ismember(Recon3.mets,mets);
    IB(IB==0) = length(mets)+1;
    samples1 = samples;
    samples2 = samples;
    samples1(samples<0 )=0; % postivitve flux
    samples2(samples>0) =0; % negative flux
    m = outmodel.S;
    %isEC = isemptycell(regexp(outmodel.subSystems,'Exchange|Transport'));
    isEC_all = Recon3.isEC;
    [EA,EB] = ismember(outmodel.rxns,Recon3.rxns);
    isEC = isEC_all(EB);

    m(:,~isEC) = 0;
    m1 = m;
    m1(m>0) = 0;%% substrate negative S
    m2 = m;
    m2(m<0) = 0; %% production positive S
    
    % consume flux  v1
    v1 = m1*samples1 + m2*samples2;
    [n_row,n_col] = size(v1);
    v1 = [v1;zeros(1,n_col)];
    v1 = v1(IB,:);
    
    % product flux v2 
    v2 = m2*samples1 + m1*samples2;
    [n_row,n_col] = size(v2);
    v2 = [v2;zeros(1,n_col)];
    v2 = v2(IB,:);
    
    % deta flux  v
    v = v1+v2;
    
    consumptionFluxMetMedian(:,i) = median(v1,2); 
    consumptionFluxMetMean(:,i) = mean(v1,2);
    consumptionFluxMetMin(:,i) = min(v1,[],2);
    consumptionFluxMetMax(:,i) = max(v1,[],2);
    productionFluxMetMedian(:,i) = median(v2,2); 
    productionFluxMetMean(:,i) = mean(v2,2);
    productionFluxMetMin(:,i) = min(v2,[],2);
    productionFluxMetMax(:,i) = max(v2,[],2);
    detaFluxMetMedian(:,i) = median(v,2); 
    detaFluxMetMean(:,i) = mean(v,2);
    detaFluxMetMin(:,i) = min(v,[],2);
    detaFluxMetMax(:,i) = max(v,[],2);
end
individual.fluxRxnsMedian = fluxRxnsMedian;
individual.fluxRxnsMean = fluxRxnsMean;
individual.fluxRxnsMin = fluxRxnsMin;
individual.fluxRxnsMax = fluxRxnsMax;
%individual.fluxFVAMin = fluxFVAMin;
%individual.fluxFVAMax = fluxFVAMax;
individual.consumptionFluxMetMedian = consumptionFluxMetMedian;
individual.consumptionFluxMetMean = consumptionFluxMetMean;
individual.consumptionFluxMetMin = consumptionFluxMetMin;
individual.consumptionFluxMetMax = consumptionFluxMetMax;

individual.productionFluxMetMedian = productionFluxMetMedian;
individual.productionFluxMetMean = productionFluxMetMean;
individual.productionFluxMetMin = productionFluxMetMin;
individual.productionFluxMetMax = productionFluxMetMax;

individual.detaFluxMetMedian = detaFluxMetMedian;
individual.detaFluxMetMean = detaFluxMetMean;
individual.detaFluxMetMin = detaFluxMetMin;
individual.detaFluxMetMax = detaFluxMetMax;
individual.modelnames = modelnames;

%% write result to ./out
if ~isdir('./out')
    mkdir('./out')
end

outheader = ['name',individual.modelnames'];
outpath = './out/';
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsMean)]],[outpath,'PQMM_fluxRxnsMean.txt']);
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsMedian)]],[outpath,'PQMM_fluxRxnsMedian.txt']);
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsMin)]],[outpath,'PQMM_fluxRxnsMin.txt']);
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsMax)]],[outpath,'PQMM_fluxRxnsMax.txt']);

writetxt([outheader;[Recon3.mets,num2cellstr(individual.consumptionFluxMetMean)]],[outpath,'PQM_consumptionFluxMetMean.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.consumptionFluxMetMedian)]],[outpath,'PQM_consumptionFluxMetMedian.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.consumptionFluxMetMin)]],[outpath,'PQM_consumptionFluxMetMin.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.consumptionFluxMetMax)]],[outpath,'PQM_consumptionFluxMetMax.txt']);

writetxt([outheader;[Recon3.mets,num2cellstr(individual.productionFluxMetMean)]],[outpath,'PQM_productionFluxMetMean.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.productionFluxMetMedian)]],[outpath,'PQM_productionFluxMetMedian.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.productionFluxMetMin)]],[outpath,'PQM_productionFluxMetMin.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.productionFluxMetMax)]],[outpath,'PQM_productionFluxMetMax.txt']);

writetxt([outheader;[Recon3.mets,num2cellstr(individual.detaFluxMetMean)]],[outpath,'PQM_detaFluxMetMean.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.detaFluxMetMedian)]],[outpath,'PQM_detaFluxMetMedian.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.detaFluxMetMin)]],[outpath,'PQM_detaFluxMetMin.txt']);
writetxt([outheader;[Recon3.mets,num2cellstr(individual.detaFluxMetMax)]],[outpath,'PQM_detaFluxMetMax.txt']);

rxnInforHeader = {'Rxns_ID','Rxns_ECnumber','Rxns_Names'};

writetxt([rxnInforHeader;[Recon3.rxns,...
    Recon3.rxnECNumbers,Recon3.rxnNames]],[outpath,'Recon3_reaction_Infor.txt']);
writetxt([{'Rxns_ID','Genes'};[Recon3.rxns,Recon3.grRules]],[outpath,'Recon3_reaction_genes.txt']);

metInforHeader = {'Mets_ID','Met_KEGG','Met_chembl','Met_PubChem',...
    'Met_HMDB','Met_name','Met_smile'};
writetxt([metInforHeader;[Recon3.mets,Recon3.metKEGGID,Recon3.metCHEBIID,...
    Recon3.metPubChemID,Recon3.metHMDBID,...
    Recon3.metNames,Recon3.metSmiles]],[outpath,'Recon3_Metabolite_Infor.txt']);    

