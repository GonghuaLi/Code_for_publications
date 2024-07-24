%% set path
addpath('./bin')
addpath('./rFASTCORMICS_for_RNA_seq_data/Scripts')
addpath('./rFASTCORMICS_for_RNA_seq_data/exampleData')
%% set overall env
pars = parse_parsfile('pars.txt');
changeCobraSolver(pars.cobrasolver)
load(pars.curatedRecon3); % a consistmodel of Recon3
Cmodel = Recon3;
%% Uptake and secration from benchmark set
vv= file2cell('./data/NCI60_flux_benchmark_Recon3(mmol L min)59.txt','\t');
exchangeRxns = vv(2:end,1);
fluxes = cell2float(vv(2:end,2:end));
vid = findRxnIDs(Recon3,exchangeRxns) >0;
exchangeRxns = exchangeRxns(vid);
fluxes = fluxes(vid,:);
addtionUptake = {'EX_his_L(e)','EX_pi(e)','EX_na1(e)','EX_fe2(e)','EX_k(e)',...
    'EX_h2o(e)','EX_Tyr_ggn(e)','EX_o2(e)'};
x1 = [-0.05;repmat(-1000,length(addtionUptake)-1,1)];
fluxes = [fluxes;repmat(x1,1,size(fluxes,2))];
exchangeRxns = [exchangeRxns; addtionUptake']; 

idx = ismember(exchangeRxns,Cmodel.rxns);
fluxes = fluxes(idx,:);
exchangeRxns = exchangeRxns(idx);


%% Read data
%data = readtable(pars.expressionFile,'ReadRowNames',1);
data = file2cell(pars.expressionFile,'\t');
rownames = data(2:end,1);%data.Properties.RowNames;
colnames = data(1,2:end);%data.Properties.VariableNames;
fpkm = cell2float(data(2:end,2:end));%(table2array(data);

%% Discretization
discretized = discretize_FPKM(fpkm, colnames);

%% disco
load('dico_rFASTCORMICS.mat')

%% setting pars
load medium_example.mat
load dico_ML.mat
already_mapped_tag = 0;
biomass_rxn = {'biomass_reaction'};
epsilon = 1e-4; 
consensus_proportion = 0.9;
unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = Cmodel.rxns(ismember(Cmodel.subSystems,unpenalizedSystems));
optional_settings.unpenalized = unpenalized;
optional_settings.func = [exchangeRxns',{'biomass_reaction','DM_atp_c_'}]; % forced additional reactions into the  model
not_medium_constrained = 'EX_tag_hs(e)';
optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium_example;% remove field if no constraint is provided

%% reconstruct model
if ~isdir('./Mat_PQM')
    mkdir('./Mat_PQM')
end
modelnames = colnames;
N = length(modelnames);
for k = 1:length(modelnames)
    [outmodel, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,k), rownames, dico, ...
        biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
    exchangeIdx = sum(abs(outmodel.S))==1 & sum(outmodel.S)==-1;
    outmodel.lb(exchangeIdx)= 0;
    for i = 1:length(exchangeRxns)
        ida = findRxnIDs(outmodel,exchangeRxns{i});
        if ida < 1
            continue;
        end
        if(fluxes(i,k) < 0)
            outmodel.lb(ida) = fluxes(i,k);
            outmodel.ub(ida) = 1000;
        else
            outmodel.lb(ida) = 0;
            outmodel.ub(ida) = 1000;
        end
    end
    outmodel.ub(findRxnIDs(outmodel,'EX_pyr(e)')) = 0;
    %FVA
    flux = FastMM_FVA_cplex(outmodel);
    outmodel.lb = flux(:,1);
    outmodel.ub = flux(:,2);
    rmRxns = outmodel.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
    outmodel = removeRxns(outmodel,rmRxns); 
    
    %save
    save(['./Mat_PQM/',modelnames{k},'.mat'],'outmodel');
    delete clone*.log;
end
