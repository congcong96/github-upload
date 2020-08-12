function badtrigfiles = intan_analysis_procedure(stimuli, probtype, flag_plot, flag_saveplot, savepath)
%multi-unit analysis for fra
if length(probtype) == 1
    multprobeopt = 0;
else
    multprobeopt = 1;
end

badtrigfiles = intan_stim_files_triggers(stimuli, multprobeopt);%multprobeopt
%return
intan_folder_raw2thresh(stimuli);
if strcmp(stimuli{1}, 'spon')
    return
end
%stimfolder = 'E:\Congcong\Documents\stimulus\Ripple_Noise';
stimfolder = 'E:\Congcong\Documents\stimulus\thalamus';
spkfile = intan_batch_process_stim(stimfolder, stimuli);

%% plots
if flag_plot == 1
    load(spkfile)
    probtypelabel = {thresh.probe};

raster_all = [];
idx = strfind(spkfile, 'fs2000');
rasterfilename = [spkfile(1:idx-1) 'fs20000-raster.mat'];
for ii = 1:length(probtype)
    %% get probinfo and deepest depth
    [probinfo] = neuronexus_prob(probtype(ii));
    probidx = probinfo.posi_idxdepth;
    probidxintan = probinfo.posi_intan;
    prob_x = probinfo.posi_x;
    idxthresh = contains(probtypelabel, probtype{ii});
    threshorig = thresh;
    thresh = threshorig(idxthresh);
    basename = thresh(1).file(1:end-9); 

    %% fra/fra10
    if  contains(stimuli,'fra')
        tcorig = tc;
        tc = tcorig(idxthresh);
        if ii == 1
            threshorig(1).sig_psth = [];
            threshorig(1).position = [];
            threshorig(1).bf = [];
        end
        raster = calculate_raster(thresh,trigger,tcparams);
        tmp = {thresh.chan};
        [raster.chan] = tmp{:};
        [~, idx] = sort([tmp{:}]);
        raster = raster(idx);
        
        figure('Renderer', 'painters', 'Position', [30 30 1000 1000]);
        for i=1:64
            f = find(i == probidx);
            chan = probidxintan(f);
            subplot(8,8,i)
            [sig(i)] = plot_psth(raster(chan+1));
            raster(chan+1).probe = probtype{ii};
            if ~isfield(thresh, 'position')|| isempty(thresh(chan+1).position)
                deepest = thresh(chan+1).depth;
                deep = deepest-probinfo.posi_depth(f);
                thresh(chan+1).position = [prob_x(f), deep];
                thresh(chan+1).sig_psth = sig(i);
            else
                deep = thresh(chan+1).position(2);
            end
            title(sprintf('A%d: %dum, %d',chan,deep,sig(i)))
        end
        if flag_saveplot
            saveas(gcf, fullfile(savepath, [basename, 'psth.png']))
        end
        close
        raster_all = [raster_all, raster];
        raster = raster_all;
        save(rasterfilename, 'raster')
        
        position = cell2mat({thresh.position}');
        [~,idx]=sortrows(position, [1 2], {'descend', 'ascend'});
        threshtmp = thresh(idx);
        
        
        figure('Renderer', 'painters', 'Position', [30 30 1000 1000]);
        if strcmp(stimuli,'fra10')
            %taxis =
        end
        for i=1:64
            chan = threshtmp(i).chan;
            deep = threshtmp(i).position(2);
            
            subplot(8,8,i)
            bf = tuning_curve_plot(tc(chan+1));
            threshtmp(i).bf = bf; %maximum of projection to frequency axis
            %title(sprintf('A%d: %dum deep, %d',chan,deep,sig(i)))
            title(sprintf('A%d: %dum %d',chan,deep, sig(i)))
        end
        if flag_saveplot
            saveas(gcf, fullfile(savepath, [basename, 'fra.png']))
        end
        close
        
        threshorig(idxthresh) = threshtmp;
        thresh = threshorig;
        tc = tcorig;
        if ii == length(probtype)
            save(spkfile, 'tcparams','thresh', '-append');
        end
    elseif flag_plot == 1 && strcmp(stimuli,'dmr')
        
        savefolder = 'E:\Congcong\Documents\emsemble_thalamus\figure\multiunit';
        savename = sprintf('%s-site%d-%dum-%s-dmr-strf', thresh(1).exp, thresh(1).site,thresh(1).depth, thresh(1).atten);
        savepath = fullfile(savefolder, savename);
        plot_strf_NH(strf(idxthresh), trigger, probinfo, thresh(1).depth, savepath)
        thresh = threshorig;
    elseif flag_plot == 1 && strcmp(stimuli,'dmrrep')
        
        savefolder = 'E:\Congcong\Documents\emsemble_thalamus\figure\multiunit';
        savename = sprintf('%s-site%d-%dum-%s-dmrrep-raster', thresh(1).exp, thresh(1).site,thresh(1).depth, thresh(1).atten);
        savepath = fullfile(savefolder, savename);
        raster_tmp = raster(idxthresh);
        plot_raster_mu(raster_tmp, probinfo, thresh(1).depth, savepath);
        thresh = threshorig;
     end
end
end