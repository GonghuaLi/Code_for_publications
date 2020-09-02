pars = parse_parsfile('pars.txt');
changeCobraSolver(pars.cobrasolver)
load(pars.curatedRecon3)
load('./data/lb_ub_for_PQMM');
vv = file2cell(pars.exchangeFile,'\t');
%%
vv= file2cell('./data/NCI60_flux_benchmark_Recon3(mmol L min)59.txt','\t');
bloodrxns = vv(2:end,1);
fluxes = cell2float(vv(2:end,2:end));
vid = findRxnIDs(Recon3,bloodrxns) >0;
bloodrxns = bloodrxns(vid);
fluxes = fluxes(vid,:);
%% uptake 
addtionUptake = {'EX_his_L(e)','EX_pi(e)','EX_na1(e)','EX_fe2(e)','EX_k(e)',...
    'EX_h2o(e)','EX_Tyr_ggn(e)','EX_o2(e)'};
x1 = [-0.05;repmat(-1000,length(addtionUptake)-1,1)];
fluxes = [fluxes;repmat(x1,1,size(fluxes,2))];
bloodrxns = [bloodrxns; addtionUptake']; 
%%
%%
N = size(lb,2);
if ~isdir('./Mat_PQM')
    mkdir('./Mat_PQM')
end

modelnames = sampleInfo(:,1);

%%
for r = 1:ceil(N/120)

outmodels = {};
%for k = 1:N
for k = (r-1)*120+1:min(N,r*120)
    outmodels{k} = Recon3;
    outmodels{k}.lb = lb(:,k);
    outmodels{k}.ub = ub(:,k);
    exchangeIdx = sum(abs(outmodels{k}.S))==1 & sum(outmodels{k}.S)==-1;
    outmodels{k}.lb(exchangeIdx)= 0;
    for i = 1:length(bloodrxns)
        ida = findRxnIDs(outmodels{k},bloodrxns{i});
        if ida < 1
            continue;
        end
        if(fluxes(i,k) < 0)
            outmodels{k}.lb(ida) = fluxes(i,k);
            outmodels{k}.ub(ida) = 1000;
        else
            outmodels{k}.lb(ida) = 0;
            outmodels{k}.ub(ida) = 1000;
        end
    end
    outmodels{k}.ub(findRxnIDs(outmodels{k},'EX_pyr(e)')) = 0;
    outmodels{k}.ub(findRxnIDs(outmodels{k},'EX_nh4(e)')) = 0;

end

%multiple threading FVA
% using parpool, you should start matlab : matlab -nodisplay  or matlab
%                do not use -nojvm option
CoreNum = pars.numfva;
if isempty(gcp('nocreate')) 
    parpool('local',CoreNum);
else  
    disp('matlab pool already started::reopen matlab pool');
    delete(gcp('nocreate'));
    parpool('local',CoreNum); 
end

%parfor k = 1:N
parfor k = (r-1)*120+1:min(N,r*120)
    %flux = FastMM_FVA(outmodels{k});
    flux = FastMM_FVA_cplex(outmodels{k});
    outmodels{k}.lb = flux(:,1);
    outmodels{k}.ub = flux(:,2);
    rmRxns = outmodels{k}.rxns(flux(:,1)>-1e-9 & flux(:,2)<1e-9);
    outmodels{k} = removeRxns(outmodels{k},rmRxns);  
end
delete(gcp('nocreate'));

% save
%for k = 1:N
for k = (r-1)*120+1:min(N,r*120)
    outmodel = outmodels{k};
    save(['./Mat_PQM/',modelnames{k},'.mat'],'outmodel');
end
eval('clear outmodels');
end



