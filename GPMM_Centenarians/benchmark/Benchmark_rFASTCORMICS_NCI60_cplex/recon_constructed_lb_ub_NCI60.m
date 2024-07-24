
%[lb,ub,outmodel] =  get_lb_ub_rxns_by_PQM_from_matrix(parsfile);
addpath('./bin');
pars = parse_parsfile('pars.txt');
changeCobraSolver(pars.cobrasolver)
load(pars.curatedRecon3)
load('./data/lb_ub_for_PQMM');
modelnames = sampleInfo(:,1);
%% uptake 
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
%%
if ~isdir('./Mat_PQM')
    mkdir('./Mat_PQM')
end

N = length(modelnames);
for k = 1:length(modelnames)
    outmodel = Recon3;
    outmodel.lb = lb(:,k);
    outmodel.ub = ub(:,k);
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
    xid = abs(outmodel.lb) < 1e-9 & abs(outmodel.ub) < 1e-9;
    outmodel = removeRxns(outmodel,outmodel.rxns(xid));
    flux = FastMM_FVA_cplex(outmodel);
    outmodel.lb = flux(:,1);
    outmodel.ub = flux(:,2);
    rmRxns = outmodel.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
    outmodel = removeRxns(outmodel,rmRxns); 
    
    %save
    save(['./Mat_PQM/',modelnames{k},'.mat'],'outmodel');
end
