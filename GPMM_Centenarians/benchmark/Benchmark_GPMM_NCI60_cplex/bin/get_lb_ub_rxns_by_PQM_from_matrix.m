function [lb,ub,outmodel] =  get_lb_ub_rxns_by_PQM_from_matrix(parsfile)
%% read parsfile
pars = parse_parsfile(parsfile);
%%
load(pars.curatedRecon3);
changeCobraSolver(pars.cobrasolver);
load('./data/mRNA_base');
%%
if strcmp(pars.expressionType,'FPKM')
    ratio_infor = file2cell('./data/ratio_FPKM.txt','\t');
    base_rna_expr = mRNA_base.fpkm_expr;
    base_rna_gene = mRNA_base.fpkm_Gene;
    base_rna_tissues = mRNA_base.fpkm_tissues;
elseif strcmp(pars.expressionType,'RMA')
    ratio_infor = file2cell('./data/ratio_RMA.txt','\t');
    base_rna_expr = mRNA_base.rma_expr;
    base_rna_gene = mRNA_base.rma_Gene;
    base_rna_tissues = mRNA_base.fpkm_tissues;
elseif strcmp(pars.expressionType,'RSEM')
    ratio_infor = file2cell('./data/ratio_RSEM.txt','\t');
    base_rna_expr = mRNA_base.rsem_expr;
    base_rna_gene = mRNA_base.rsem_Gene;
    base_rna_tissues = mRNA_base.rsem_tissues;
else
    error('PQMM just supports three type of expression data: RMA,FPKM and RSEM.\n');
end
ratio = cell2float(ratio_infor(:,2));
geneNames = ratio_infor(:,1);
expr_infor = file2cell(pars.expressionFile,'\t');
%%
sampleInfo = expr_infor(1,2:end)';
%%
expr = cell2float(expr_infor(2:end,2:end));
expr_geneNames = expr_infor(2:end,1);
% normalize expr base
if strcmp(pars.standardization,'yes')
    [IA,IB] = ismember(expr_geneNames,base_rna_gene);
    IB  = IB(IA);
    common_base_expr = base_rna_expr(IB,:);
    %if ~strcmp(pars.expressionType,'RMA')
        
    basetmp = log2(common_base_expr);
    basetmp(basetmp < -13.29) = nan;
    log2_mu = median(nanmean(basetmp));
    log2_std = median(nanstd(basetmp));
        
    log2_expr = log2(expr);
    log2_expr(log2_expr < -13.29) = nan;
    log2_expr_comm = log2_expr(IA,:);
    log2_expr_mu =  nanmean(log2_expr_comm);
    log2_expr_std =  nanstd(log2_expr_comm);
        
    log2_expr = ((log2_expr - repmat(log2_expr_mu,size(expr,1),1))./repmat(log2_expr_std,size(expr,1),1))...
        *log2_std + log2_mu;
    expr = 2.^log2_expr;
    expr(isnan(expr)) = 0;
end      
N = size(expr,2);
%%
[IA,IB] = ismember(geneNames,expr_geneNames);
geneNames = geneNames(IA);
ratio = ratio(IA);
IB = IB(IA);
expr = expr(IB,:);
if strcmp(pars.expressionType,'RMA')
    expr = 2.^expr;
end
predExpr =expr.*repmat(ratio,1,N);
%% 
lb = zeros(length(Recon3.rxns),N);
ub = zeros(length(Recon3.rxns),N);
%Recon3.Kcat(findRxnIDs(Recon3,'PDHm'),:) = 2;  % recurated PDHm
%%
if strcmp(pars.uptakeUnite,'uM/100ml/min')
    unitfactor = 0.006;
elseif strcmp(pars.uptakeUnite,'mM/L/min')
    unitfactor = 0.00006;
else
    error('PQMM just supports two type uptake units: uM/100ml/min or mM/L/min.\n');
end
for i = 1:N
    enzyme_expr = recon_get_enzyme_expr(geneNames,predExpr(:,i),Recon3.geneSymbol,0);
    [outmodel,par] = recon_quantitative_model_nofva(Recon3,enzyme_expr,unitfactor );
    lb(:,i) = outmodel.lb;
    ub(:,i) = outmodel.ub;
end
save('./data/lb_ub_for_PQMM','lb','ub','sampleInfo');
    
function enzyme_expr = recon_get_enzyme_expr(genes,expr,enzyme_genes,cutoff)
idx = expr<cutoff;
tissue_expr = expr;
tissue_expr(idx) = 0;
[ia,ib] = ismember(enzyme_genes,genes);
enzyme_expr = zeros(size(ia));
for i = 1:length(ib)
    if ib(i) ==0
       enzyme_expr(i) =  -1000;
    else
       enzyme_expr(i) = tissue_expr(ib(i));
    end
end