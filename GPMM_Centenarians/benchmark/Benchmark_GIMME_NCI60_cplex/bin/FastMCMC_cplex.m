function outmodel = FastMCMC_cplex(inmodel,numPoints,numSteps,outfile,varargin)

tic;
t = clock;
ud =0;
% mat model
lb = inmodel.lb;
ub = inmodel.ub;
b = inmodel.b;
c = zeros(length(inmodel.rxns),1);
[bed,ind,val] = sparse_to_csr(inmodel.S);
bed = bed-1;
ind = ind-1;
a1 = ceil(rand(1,1)*10000);
a2 = ceil(rand(1,1)*10000);
b1 = ceil(rand(1,1)*10000);
b2 = ceil(rand(1,1)*10000);
inmatname = ['modelIn_cplex_',num2str(a1),num2str(a2),num2str(ceil(t(6)*100)),'.mat'];
save(inmatname,'lb','ub','b','c','bed','ind','val');


% null_S
nullS_varname = ['nullS',num2str(b1),num2str(b2),num2str(ceil(t(6)*100))];
tmp_nullS = null(full(inmodel.S));
eval([nullS_varname,' =  tmp_nullS ; clear tmp_nullS'])
save(nullS_varname,nullS_varname);
nullSname = [nullS_varname,'.mat'];  

if length(varargin)>2
    ud =1;
    mTol = varargin{1};
    uTol = varargin{2};
    dTol = varargin{3};
end

if length(outfile)>4 && strcmpi(outfile(end-3:end),'.mat')
    cout = outfile;
else
    cout = [outfile,'.mat'];
end

% rum mcmc
if ud ==0
    mcmc_cmd = ['mcmc_cplex_mkl_prepared_null ',inmatname,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout]
else
    mcmc_cmd = ['mcmc_cplex_mkl_prepared_null ',inmatname,' ',nullSname,' ',num2str(numPoints),...
        ' ',num2str(numSteps),' ',cout,' ',mTol,' ',uTol,' ',dTol]
end
system(mcmc_cmd);

load(cout);
outmodel = inmodel;
outmodel.points = points;

delete(nullSname);
delete(inmatname);
toc






