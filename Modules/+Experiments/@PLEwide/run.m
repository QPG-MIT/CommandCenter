function run( obj,status,managers,ax )
    % Main run method (callback for CC run button)
    obj.abort_request = false;
    status.String = 'Experiment started';
    drawnow;
    
    rangeWL = eval(obj.ScanRange);
    rangeFreq = Sources.TunableLaser_invisible.c./rangeWL; % change wavelength (nm) to frequency (THz)
    obj.PM.set_average_count(1);
    PWavg_num = 20;
    PWall = zeros(length(rangeWL),PWavg_num);
    
    switch obj.Detection_Type
        case 'Spectrometer'
            panel = ax.Parent; delete(ax); % For spectrometer detection, we need two axes for plotting
            ax(1) = subplot(1,2,1,'parent',panel);
            ax(2) = subplot(1,2,2,'parent',panel);
            obj.takeSpec.run(status,managers,ax(2)) % Take Dummy Spectra to get spectrometer pixel number
            spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);
            scanH = imagesc(ax(1),spectrum.wavelength,rangeWL,NaN(length(rangeWL),length(spectrum.wavelength))); % ,'parent',managers.handles.axImage
            set(ax(1),'ydir','normal');
        case 'APD'
            scanH = plot(ax,rangeWL,NaN(1,length(rangeWL)),'LineWidth',1.5);
            xlabel(ax,'Excitation Wavelength (nm)')
            ylabel(ax,'APD counts')
    end
    
    for i = 1:length(rangeWL)
        assert(~obj.abort_request,'User aborted.');
        
%         obj.resLaser.TuneCoarse(rangeFreq(i)); % Tune solstis or EMM
%         obj.data.laser_wavelength(i) = Sources.TunableLaser_invisible.c/obj.resLaser.getFrequency; % Get exact wavelength from wavemeter


        obj.PM.set_wavelength(round(rangeWL(i))); % set powermeter wavelength
        for j =1:PWavg_num % Take power data
            try
                PWall(i,j) = obj.PM.get_power('MW');
            catch err
                warning(err.message)
                PWall(i,j) = NaN;
                continue
            end
        end
        obj.data.laser_power(i) = nanmean(PWall(i,:)); % Save average power measured
        
        switch obj.Detection_Type
            case 'Spectrometer'
                obj.takeSpec.run(status,managers,ax(2)) % Run spectrum
                spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);  % Get spectrum data
                obj.data.spec_wavelength = spectrum.wavelength;
                scanH.CData(i,:) = spectrum.intensity;
                title(ax(2),sprintf('Spectra %i of %i',i,length(rangeWL)));
                drawnow
            case 'APD'
                APD_counts = obj.counter.singleShot(obj.APD_dwell)*obj.APD_dwell/1000;
                scanH.YData(i) = APD_counts;
%                 plot(ax,rangeWL(i),obj.data.APD_counts(i))
                drawnow
        end

    end
    
    switch obj.Detection_Type
        case 'Spectrometer'
            obj.data.spec_intensity = scanH.CData;
        case 'APD'
            obj.data.APD_counts = scanH.YData; 
            obj.data.APD_dwelltime = obj.APD_dwell; 
    end
    obj.data.set_wavelength = rangeWL;
    
    % Edit this to include meta data for this experimental run (saved in obj.GetData)
    obj.meta.prefs = obj.prefs2struct;
    obj.meta.position = managers.Stages.position; % Save current stage position (x,y,z);
    
    
    try
        % EXPERIMENT CODE %
    catch err
    end
    % CLEAN UP CODE %
    if exist('err','var')
        % HANDLE ERROR CODE %
        rethrow(err)
    end
end
