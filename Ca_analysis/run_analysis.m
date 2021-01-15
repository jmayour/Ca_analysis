function run_analysis(mainfile,parameterfile)

parameters = readtable(parameterfile,'Range','B:B'); % call parameters
[~, savename, ~] = fileparts(mainfile); % naming of output excel file

% run analysis
calcium_analysis(mainfile,savename,parameters{2,1},parameters{3,1},parameters{5,1},...
    parameters{6,1},parameters{7,1}./100,parameters{10,1},parameters{11,1},parameters{12,1},...
    parameters{13,1},parameters{8,1},parameters{14,1},parameters{15,1},parameters{9,1},...
    parameters{1,1},parameters{4,1});

end
