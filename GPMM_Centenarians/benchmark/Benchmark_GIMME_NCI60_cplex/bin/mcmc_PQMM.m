%% run individual mcmc
%normal
addpath('./bin');
load('./data/lb_ub_for_PQMM');
modelnames = sampleInfo(:,1);
pars = parse_parsfile('pars.txt');
changeCobraSolver(pars.cobrasolver);
for i = 1:length(modelnames)
     tname = ['./Mat_PQM/MCMC',modelnames{i},'_1.mat'];
     if exist(tname)
       continue;
     end
    load(['./Mat_PQM/',modelnames{i},'.mat']);
    idatp = findRxnIDs(outmodel,'DM_atp_c_');
    outmodel.lb(idatp) = outmodel.ub(idatp)*0.5;%just for GIMME 0.5
    idbiomass  = findRxnIDs(outmodel,'biomass_reaction');
    sol= optimizeCbModel(outmodel);
    outmodel.lb(idbiomass)  = sol.f*0.5;%just for GIMME 0.5
     outmodel.ub(idbiomass)  = sol.f;
    warmupPts= createHRWarmup(outmodel);
     if size(warmupPts,1) < 2
        continue;
     end
    ACHRSampler(outmodel,warmupPts,['./Mat_PQM/MCMC',modelnames{i}],1,pars.numPoints,1000);
end
exit;