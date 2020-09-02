function [score0,out,scores] = check_aa(line,ExprStat)
out = [];

        line = regexprep(line,'x\(','');
        line = regexprep(line,'(\d)\)','$1');
        line = regexprep(line,' ','');
        score0 = -1000;
        scores = -1000*ones(100,1);
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
                    score0 = max(score0,ExprStat(str2num(line(k1:k2))));
                    out = [out;i,score0];
                elseif i == length(line) & opr == '&'
                    score0 = min(score0,ExprStat(str2num(line(k1:k2))));
                    out = [out;i,score0];
                else
                    continue;
                end
            elseif line(i) == '&' | line(i) == '|'
                if line(i+1) =='('
                    opr2 = line(i);
                    idscore = idscore+1
                    continue;
                end
                if opr == '|'
                    scores(idscore) = max(scores(idscore),ExprStat(str2num(line(k1:k2))));
                    if opr2 == ' '
                        score0 = max(score0,scores(idscore));
                        out = [out;i,score0];
                    end
                else
                    scores(idscore) = min(scores(idscore),ExprStat(str2num(line(k1:k2))));
                    if opr2 == ' '
                        score0 = min(score0,scores(idscore));
                        out = [out;i,score0];
                    end
                end
                opr = line(i);
                k1 = i+1;
                k2 = i+1;
            elseif line(i) == ')'                  
                if opr == '|'
                    scores(idscore) = max(scores(idscore),ExprStat(str2num(line(k1:k2))));
                else
                    scores(idscore) = min(scores(idscore),ExprStat(str2num(line(k1:k2))));
                end
                if opr2 == '&'
                    score0 = min(score0,scores(idscore));
                else
                    score0 = max(score0,scores(idscore));
                end
                opr2 = ' ';
                out = [out;i,score0];
            end
            
        end
    end
