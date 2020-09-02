function [rxnState,core_rxns] = get_rxns_expr_byrule(modelin,enzyme_expr)
% reconstruction Cb model through a given expressed genelist.
% by Gonghua Li
% 2014-04-17
%modelin = liver;
%genstat = gene_stat;
%genelist = modelin.genes(genstat>0);
genes = modelin.genes;
%genesEntrez = regexp(genes,'([^\.]+)\.\d','tokens');
%genesEntrez = [genesEntrez{:}];
%genesEntrez = [genesEntrez{:}];
rules = modelin.rules;
%%ExprState = ismember(genesEntrez,genelist);
ExprStat = enzyme_expr;%%ismember(genes,genelist);
rxnState  = zeros(length(modelin.rxns),1);
%%
for j = 1:length(rules)
    if isempty(rules{j})
        rxnState(j) = -1000;
    else
        line = rules{j};
        line = regexprep(line,'x\(','');
        line = regexprep(line,'(\d)\)','$1');
        line = regexprep(line,' ','');
        score0 = 0;
        scores = zeros(1000,1);
        k1 = 1;
        k2 = 1;
        idscore = 1;
        opr = '|';
        opr2 = ' ';
        for i = 1:length(line)
            if line(i) == '('
                k1 = i+1;
                k2 = i+1;
                opr = '|';
            elseif line(i) >= '0' & line(i) <= '9'
                k2 = i;
                if i == length(line) & opr == '|';
                    score0 = sum2(score0,ExprStat(str2num(line(k1:k2))));
                elseif i == length(line) & opr == '&'
                    score0 = min(score0,ExprStat(str2num(line(k1:k2))));
                else
                    continue;
                end
            elseif line(i) == '&' | line(i) == '|'
                if line(i+1) =='('
                    opr2 = line(i);
                    idscore = idscore+1;
                    continue;
                end
                if opr == '|'
                    scores(idscore) = sum2(scores(idscore),ExprStat(str2num(line(k1:k2))));
                    if opr2 == ' '
                        score0 = sum2(score0,scores(idscore));
                    end
                else
                    scores(idscore) = min(scores(idscore),ExprStat(str2num(line(k1:k2))));
                    if opr2 == ' '
                        score0 = min(score0,scores(idscore));
                    end
                end
                opr = line(i);
                k1 = i+1;
                k2 = i+1;
            elseif line(i) == ')'                  
                if opr == '|'
                    scores(idscore) = sum2(scores(idscore),ExprStat(str2num(line(k1:k2))));
                else
                    scores(idscore) = min(scores(idscore),ExprStat(str2num(line(k1:k2))));
                end
                if opr2 == '&'
                    score0 = min(score0,scores(idscore));
                else
                    score0 = sum2(score0,scores(idscore));
                end
                opr2 = ' ';
            end  
        end
        rxnState(j) = score0;
        if rxnState(j) == 0
            rxnState(j) = -1000;
        end
    end
end
      
core_rxns = find(rxnState > 1e-6);
%%rxnRemoves = modelin.rxns(rxnState<0.01);
%%modelout = changeRxnBounds(modelin,rxnRemoves,0,'b');
%%modelout.geneStates = ExprState;