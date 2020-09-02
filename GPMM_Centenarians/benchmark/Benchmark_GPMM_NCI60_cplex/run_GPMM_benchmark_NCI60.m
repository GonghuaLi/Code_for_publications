addpath('./bin');
Reduce_recon3_benchmark_NCI60;
get_lb_ub_rxns_by_PQM_from_matrix('pars.txt');
recon_PQMM_multi_benchmark_NCI60;
multi_threading_mcmc;



