datafolder = 'E:\Congcong\Documents\emsemble_thalamus\2020-01-27-CH';

%% FRA
stimuli = {'fra'};
probtype = {'H31x64'};
flag_plot = 1;
flag_saveplot = 1;
savepath = datafolder;
 
frafolders = dir(fullfile(datafolder,'*_fra_*'));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath);%
end

%% FRA10
stimuli = {'fra10'};
probtype = {'H31x64'};
flag_plot = 1;
flag_saveplot = 1;
savepath = datafolder;
 
frafolders = dir(fullfile(datafolder,'*_fra10_*'));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath);%
end

%% dmr
stimuli = {'dmr'};
probtype = {'H31x64'};
flag_plot = 1;
flag_saveplot = 1;
savepath = datafolder;
 
frafolders = dir(fullfile(datafolder,'*_dmr_*'));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath);%
end

%% dmrrep
stimuli = {'dmrrep'};
probtype = {'H31x64'};
flag_plot = 1;
flag_saveplot = 1;
savepath = datafolder;
 
frafolders = dir(fullfile(datafolder,'*_dmrrep_*'));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath);%
end

%% spon
stimuli = {'spon'};
probtype = {'H31x64'};
flag_plot = 1;
flag_saveplot = 1;
savepath = datafolder;
 
frafolders = dir(fullfile(datafolder,'*_spon_*'));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath);%
end
%% plot bf v.s. recording depth
frafolders = dir(fullfile(datafolder,'*_fra_*'));
depth = cell(1,length(frafolders));
bf = cell(1,length(frafolders));
for i = 1:length(frafolders)
    cd(fullfile(datafolder, frafolders(i).name))
    framat = dir('*fra.mat');
    load(framat.name)
    for j = 1:64
        if thresh(j).sig_psth == 0
            continue
        end
        depth{i} = [depth{i} thresh(j).position];
        bf{i} = [bf{i} thresh(j).bf];
    end
    scatter(depth{i}, bf{i}/1000,'filled')
    hold on
end

set(gca, 'YScale', 'log')
xlabel('depth/um')
ylabel('bf/kHz')
yticks([2, 4, 8, 16, 32])
yticklabels({'2', '4', '8', '16', '32'})