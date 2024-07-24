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
        outmodels{k}.lb(findRxnIDs(outmodels{k},bloodrxns{i})) = exlb(i);
        outmodels{k}.ub(findRxnIDs(outmodels{k},bloodrxns{i})) = exub(i);
    end
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



