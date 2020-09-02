% summary PQMM results
pars = parse_parsfile('pars.txt');
% get_individual_result V0.1
load(pars.curatedRecon3);
load('./data/lb_ub_for_PQMM.mat');
%changeCobraSolver('gurobi5');
changeCobraSolver('cplex');
[num_mets,num_rxns] = size(Recon3.S);

modelnames = sampleInfo(:,1);
N_samples = length(modelnames);

fluxRxnsMedian = zeros(num_rxns,N_samples);
fluxRxnsMean = zeros(num_rxns,N_samples);
fluxRxnsMin = zeros(num_rxns,N_samples);
fluxRxnsMax = zeros(num_rxns,N_samples);
fluxFVAMin = zeros(num_rxns,N_samples);
fluxFVAMax = zeros(num_rxns,N_samples);
fluxRxnsCV = zeros(num_rxns,N_samples);
fluxRxnsEntropy = zeros(num_rxns,N_samples);
pars = parse_parsfile('pars.txt');
for i = 1:length(modelnames)
    i
    if exist(['./Mat_PQM/MCMC',modelnames{i},'_1.mat'],'file')   
        samples = loadSamples(['./Mat_PQM/MCMC',modelnames{i}],1,pars.numPoints);
    else
         continue;
    end
    load(['./Mat_PQM/',modelnames{i},'.mat']);
    rxns = outmodel.rxns;
    [IA,IB] = ismember(Recon3.rxns,rxns);
    IB(IB==0) = length(rxns)+1;
    %samples = loadSamples(['./Mat_PQM/MCMC',modelnames{i}],1,10000);
    % rm all zero points%
    samples = samples(:,sum(abs(samples),1) > 1e-5);
    %
    medianFlux = median(samples,2);
    meanFlux = mean(samples,2);
    minFluxMCMC = min(samples,[],2);
    maxFluxMCMC = max(samples,[],2);
    tcv = cal_CV(samples);
    %tentropy = calc_entropy(samples,3.23,200);
     tentropy = calc_entropy2_1(samples,200);
    medianFlux = [medianFlux;0];
    meanFlux = [meanFlux;0];
    minFluxMCMC = [minFluxMCMC;0];
    maxFluxMCMC = [maxFluxMCMC;0];
    %fluxMinMax = [fluxMinMax;0,0];
    tcv = [tcv;0];
    tentropy = [tentropy;0];
    
    fluxRxnsMedian(:,i) = medianFlux(IB);
    fluxRxnsMean(:,i) = meanFlux(IB);
    fluxRxnsMin(:,i) = minFluxMCMC(IB);
    fluxRxnsMax(:,i) = maxFluxMCMC(IB);
    %fluxFVAMin(:,i) = fluxMinMax(IB,1);
    %fluxFVAMax(:,i) = fluxMinMax(IB,2);
    fluxRxnsCV(:,i) = tcv(IB);
    fluxRxnsEntropy(:,i) = tentropy(IB);
end
individual.fluxRxnsMedian = fluxRxnsMedian;
individual.fluxRxnsMean = fluxRxnsMean;
individual.fluxRxnsMin = fluxRxnsMin;
individual.fluxRxnsMax = fluxRxnsMax;
%individual.fluxFVAMin = fluxFVAMin;
%individual.fluxFVAMax = fluxFVAMax;
individual.fluxRxnsCV = fluxRxnsCV;
individual.fluxRxnsEntropy = fluxRxnsEntropy;

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
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsCV)]],[outpath,'PQMM_fluxRxnsCV.txt']);
writetxt([outheader;[Recon3.rxns,num2cellstr(individual.fluxRxnsEntropy)]],[outpath,'PQMM_fluxRxnsEntropy.txt']);

rxnInforHeader = {'Rxns_ID','Rxns_ECnumber','Rxns_Names'};

writetxt([rxnInforHeader;[Recon3.rxns,...
    Recon3.rxnECNumbers,Recon3.rxnNames]],[outpath,'Recon3_reaction_Infor.txt']);
writetxt([{'Rxns_ID','Genes'};[Recon3.rxns,Recon3.grRules]],[outpath,'Recon3_reaction_genes.txt']);

metInforHeader = {'Mets_ID','Met_KEGG','Met_chembl','Met_PubChem',...
    'Met_HMDB','Met_name','Met_smile'};
writetxt([metInforHeader;[Recon3.mets,Recon3.metKEGGID,Recon3.metCHEBIID,...
    Recon3.metPubChemID,Recon3.metHMDBID,...
    Recon3.metNames,Recon3.metSmiles]],[outpath,'Recon3_Metabolite_Infor.txt']);    

