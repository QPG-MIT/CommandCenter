function run( obj,status,managers,ax )
    % Main run method (callback for CC run button)
    obj.abort_request = false;
    status.String = 'Experiment started';
    drawnow;
    
    try
        % Home rotation mount
        if ~obj.rot.Homed
            status.String = 'Homing motor'; drawnow;
            
            obj.rot.home();
            pause4Move(obj.rot, obj.motor_home_time);
        end
        
        obj.resLaser.off
        obj.data = [];
        rangeWL = eval(obj.ScanRange);
        rangeFreq = Sources.TunableLaser_invisible.c./rangeWL; % change wavelength (nm) to frequency (THz)
        rangePOL = eval(obj.PolarizationAngleRange);
        Nangles = length(rangePOL);
        obj.PM.set_average_count(1);
        PWavg_num = 25;
        PWall = zeros(length(eval(obj.ScanRange)),PWavg_num);
        
        switch obj.Detection_Type
            case 'Spectrometer'
                panel = ax.Parent; delete(ax); % For spectrometer detection, we need two axes for plotting
                ax(1) = subplot(1,2,1,'parent',panel);
                ax(2) = subplot(1,2,2,'parent',panel);
                temp_exposure = obj.takeSpec.exposure;
                obj.takeSpec.exposure = 0.001;
                obj.takeSpec.run(status,managers,ax(2)) % Take Dummy Spectra to get spectrometer pixel number
                obj.takeSpec.exposure = temp_exposure;
                spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);
                scanH = imagesc(ax(1),spectrum.wavelength,rangeWL,NaN(length(rangeWL),length(spectrum.wavelength))); % ,'parent',managers.handles.axImage
                set(ax(1),'ydir','normal');
            case 'APD'
                scanH = plot(ax,rangeWL,NaN(1,length(rangeWL)),'LineWidth',1.5);
                xlabel(ax,'Excitation Wavelength (nm)')
                ylabel(ax,'APD counts')
                obj.data.APD_dwelltime = obj.APD_dwell;
        end
        
        % tune wavelength
        error_count = 0;
        for i = 1:length(rangeWL)
            assert(~obj.abort_request,'User aborted.');
            try % Try wavelength tuning once more if it fails
                obj.resLaser.TuneCoarse(rangeFreq(i)); % Tune solstis or EMM
            catch err
                if strcmp(err.msgID, 'HWSERVER:empty') || strcmp(err.message,'Tuning failed')
                    try % Move to next target wavelength if it fails once again
                        obj.resLaser.TuneCoarse(rangeFreq(i)); % Tune solstis or EMM
                    catch err
                        if strcmp(err.msgID, 'HWSERVER:empty') || strcmp(err.message,'Tuning failed')
                            error_count = error_count+1;
                            obj.meta.error_count(error_count) = i;
                            obj.meta.error_msg(error_count) = err.message;
                            continue
                        else
                            rethrow(err);
                        end
                    end
                else
                    retrhow(err);
                end
            end
            pause(0.2);
            obj.resLaser.on
            pause(0.1);
            
            
            % tune laser power to desired value
            tuneLaserPower(obj,power)
            
            % sweep through angles
            for k = 1:Nangles
                
                theta = obj.angle_list(k);
                status.String = sprintf( 'Navigating to %g (%i/%i)', theta, k, Nangles); drawnow;
                obj.rot.move(theta);
                pause4Move(obj.rot, obj.motor_move_time);
                status.String = sprintf( 'Measuring at %g (%i/%i)', theta, k, Nangles); drawnow;
           
                switch obj.Detection_Type
                    case 'Spectrometer'
                        obj.takeSpec.run(status,managers,ax(2)) % Run spectrum
                        spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);  % Get spectrum data
                        obj.data.spec_wavelength = spectrum.wavelength; %
                        obj.data.spec_intensity(i,k,:) = spectrum.intensity;
                        scanH.CData(i,:) = spectrum.intensity;
                        title(ax(2),sprintf('Spectra %i of %i',(i-1)*range(POL)+k,length(rangeWL)*length(rangePOL));
                        drawnow
                    case 'APD'
                        APD_counts = obj.counter.singleShot(obj.APD_dwell)*obj.APD_dwell/1000;
                        scanH.YData(i) = APD_counts;
                        obj.data.APD_counts(i,k) = APD_counts;
                        drawnow
                end
            end
            
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
            obj.data.laser_power_all(i,:) = PWall(i,:); % Save average power measured
            obj.data.laser_wavelength(i) = Sources.TunableLaser_invisible.c/obj.resLaser.getFrequency; % Get exact wavelength from wavemeter
            obj.resLaser.off
            
        end
        obj.data.laser_setwavelength = rangeWL;
        obj.data.Detection_Type = obj.Detection_Type;
        
        % Edit this to include meta data for this experimental run (saved in obj.GetData)
        obj.meta.prefs = obj.prefs2struct;
        obj.meta.position = managers.Stages.position; % Save current stage position (x,y,z);
    
        % EXPERIMENT CODE %
    catch err
        obj.resLaser.off
    end
    % CLEAN UP CODE %
    if exist('err','var')
        % HANDLE ERROR CODE %
        rethrow(err)
    end
end


% Wait until motor stops moving, or timeout
function pause4Move(rot,maxTime)
    t = tic;
    while (rot.Moving || ~rot.Homed) && (toc(t) < maxTime)
        drawnow
        if toc(t) > maxTime
            error('Motor timed out while moving')
        end
    end
end

% tune laser power
function tuneLaserPower(obj,power)
    PWavg_num = 25;
    deltaAngle = 5;
    
    initPower = get_mean_power(obj, PWavg_num);
    initPos = obj.rotPower.position;
    power = initPower;
    steps = 1;
    
    while ( abs(Power-obj.LaserPowerTarget) > obj.LaserPowerTolerance)
        deltaPower = power-obj.LaserPowerTarget;
        deltaAngle_old = deltaAngle;
        if deltaPower < 0 % sign of power difference determines move direction
            deltaAngle = abs(deltaAngle) * -1;
        end
        if deltaAngle*deltaAngle_old < 0: % if direction changes the stepsize is halfed
            deltaAngle = deltaAngle * 0.5;
        end
        
        obj.rotPower.move_by(deltaAngle);
        pause4Move(obj.rotPower, obj.motor_move_time);
        
        power = get_mean_power(obj, PWavg_num);
        steps = steps + 1;
        
        if steps > 10
            error('Laser power tuning failed')
            obj.rotPower.move(initPos);
            pause4Move(obj.rotPower, obj.motor_move_time);
            break
        end
        
    end

end


function get_mean_power(obj,N)
    for i = 1:N
        try
            power(i) = obj.PM.get_power('MW');
        catch err
            warning(err.message)
            power(i) = NaN;
            continue
        end
    end
    return nanmean(power);
end