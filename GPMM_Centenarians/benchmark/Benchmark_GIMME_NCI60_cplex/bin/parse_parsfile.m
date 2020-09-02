function pars = parse_parsfile(parsfile)
parsInfor = file2cell(parsfile,'\n');
for i  = 1:length(parsInfor)
    thisline = parsInfor{i};
    if strcmp(thisline(1),'#') || length(thisline)<5
        continue;
    end
    m = regexp(thisline,'''([^'']+)''','tokens');
    m = m{1};
    m = m{1};
    if regexp(thisline,'^curatedRecon3')
        pars.curatedRecon3 = m;
    elseif regexp(thisline,'^cobrasolver')
        pars.cobrasolver = m;
    elseif regexp(thisline,'^isFastMM')
        pars.isFastMM = m;
    elseif regexp(thisline,'^exchangeFile')
        pars.exchangeFile = m;
    elseif regexp(thisline,'^uptakeUnite')
        pars.uptakeUnite = m;
    elseif regexp(thisline,'^expressionFile')
        pars.expressionFile = m;
    elseif regexp(thisline,'^standardization')
        pars.standardization = m;
    elseif regexp(thisline,'^expressionType')
        pars.expressionType = m;
    elseif regexp(thisline,'^sampleInfor')
        pars.sampleInfor = m;
    elseif regexp(thisline,'^degradation')
        pars.degradation = m;
    elseif regexp(thisline,'^fva_threading')
        pars.numfva = str2num(m);
    elseif regexp(thisline,'^cpu')
        pars.numCPU = str2num(m);
    elseif regexp(thisline,'^numPoints')
        pars.numPoints = str2num(m);
    end
end
        
        