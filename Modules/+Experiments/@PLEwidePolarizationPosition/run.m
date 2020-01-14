function run( obj,status,managers,ax )
    % Main run method (callback for CC run button)
    obj.abort_request = false;
    status.String = 'Experiment started';
    drawnow;
    
    stage = managers.Imaging.active_module.uses_stage;
    assert(logical(managers.Stages.check_module_str(stage)),'Stage associated with active imager is not loaded.');
    managers.Stages.setActiveModule(stage); % Make the stage active
    startingPos = managers.Stages.position;
    zpos = startingPos(3);
    
    x = eval(obj.Xcoords);
    y = eval(obj.Ycoords);
    assert(length(x) == length(y), 'x and y coord length should be the same.');
        

    try
%         % Home rotation mount
%         obj.rot = Drivers.APTMotor.instance(27251915, [0 360]);
%         obj.rotPower = Drivers.APTMotor.instance(83847313, [0 360]);
        if ~obj.rot.Homed
            status.String = 'Homing motor'; drawnow;
            obj.rot.home();
            pause4Move(obj.rot, obj.motor_home_time);
            status.String = 'motor homed'; drawnow;
        end
        if ~obj.rotPower.Homed
            status.String = 'Homing Power motor'; drawnow;
            obj.rotPower.disable
            pause(7);
            obj.rotPower.enable
            pause(1);
            obj.rotPower.home();
            pause4Homing(obj.rotPower, obj.motor_home_time);
                
            status.String = 'Power motor homed'; drawnow;
        end
        startpos = obj.rotPower.Position;
        obj.resLaser.off
        obj.data = [];
        rangeWL = eval(obj.ScanRange);
        rangeFreq = Sources.TunableLaser_invisible.c./rangeWL; % change wavelength (nm) to frequency (THz)
        rangePOL = eval(obj.PolarizationAngleRange);
        Nangles = length(rangePOL);
        obj.PM.set_average_count(1);
        PWavg_num = 25;
        Nsteps = 26;
        PWall = zeros(length(eval(obj.ScanRange)),PWavg_num);
        
        switch obj.Detection_Type
            case 'Spectrometer'
                panel = ax.Parent; delete(ax); % For spectrometer detection, we need two axes for plotting
                ax(1) = subplot(2,2,1,'parent',panel);
                ax(2) = subplot(2,2,2,'parent',panel);
                temp_exposure = obj.takeSpec.exposure;
                obj.takeSpec.exposure = 0.001;
                obj.takeSpec.run(status,managers,ax(2)) % Take Dummy Spectra to get spectrometer pixel number
                obj.takeSpec.exposure = temp_exposure;
                spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);
                scanH = imagesc(ax(1),spectrum.wavelength,rangeWL,NaN(length(rangeWL),length(spectrum.wavelength))); % ,'parent',managers.handles.axImage
                set(ax(1),'ydir','normal');
            case 'APD'
                panel = ax.Parent; delete(ax); % For spectrometer detection, we need two axes for plotting
                tempax = subplot(2,2,1,'parent',panel);
                ax(2) = subplot(2,2,2,'parent',panel);
                ax(1) = polaraxes(panel,'Units',tempax.Units,'Position',tempax.Position);
                delete(tempax);
                scanH = plot(ax(2),rangeWL,NaN(1,length(rangeWL)),'LineWidth',1.5,'MarkerSize',5);
                xlabel(ax(2),'Excitation Wavelength (nm)')
                ylabel(ax(2),'APD counts')
                scanP = polarplot(ax(1),deg2rad(rangePOL)*2,NaN(1,length(rangePOL)),'LineWidth',1.5,'MarkerSize',5);
%                 xlabel(ax(1),'Polarization Angle (deg)')
%                 ylabel(ax(1),'APD counts')
                obj.data.APD_dwelltime = obj.APD_dwell;
        end
                ax(3) = subplot(2,2,3,'parent',panel);
                scanPW = plot(ax(3),rangeWL,NaN(1,length(rangeWL)),'LineWidth',1.5,'MarkerSize',5);
                xlabel(ax(3),'Excitation Wavelength (nm)')
                ylabel(ax(3),'Tuned Laser Power (uW)')
                
                ax(4) = subplot(2,2,4,'parent',panel);
                scanPT = plot(ax(4),NaN(1,Nsteps),NaN(1,Nsteps),'LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor','b');
                xlabel(ax(4),'Position (a.u.)')
                ylabel(ax(4),'Tuned Laser Power (uW)')
        
        % tune wavelength
%         error_count = 0;
        for i = 1:length(rangeWL)
            assert(~obj.abort_request,'User aborted.');
            
            
            current_freq = rangeFreq(i);
            error_flag = 5;
            while error_flag
                assert(~obj.abort_request,'User aborted.');
                try
                    obj.resLaser.TuneCoarse(current_freq); % Tune solstis or EMM
                    break
                catch err
                    current_freq = current_freq-0.01;
                    error_flag = error_flag-1;
                    msg = sprintf('Tuning error flag value: %i',error_flag)
                    pause(20);
                end
            end
            if ~error_flag
                msg = 'Tuning failed. Continue..'
                continue    
            end
           
            
            
            
            
%             try % Try wavelength tuning once more if it fails
%                 obj.resLaser.TuneCoarse(rangeFreq(i)); % Tune solstis or EMM
%             catch err
% %                 if strcmp(err.msgID, 'HWSERVER:empty') || strcmp(err.message,'Tuning failed')
%                     try % Move to next target wavelength if it fails once again
%                         obj.resLaser.TuneCoarse(rangeFreq(i)-0.01); % Tune solstis or EMM
%                     catch err
% %                         if strcmp(err.msgID, 'HWSERVER:empty') || strcmp(err.message,'Tuning failed')
% %                             error_count = error_count+1;
% %                             obj.meta.error_count(error_count) = i;
% %                             obj.meta.error_msg(error_count) = err.message;
%                             continue
% %                         else
% %                             rethrow(err);
% %                         end
%                     end
% %                 else
% %                     retrhow(err);
% %                 end
%             end
            pause(0.2);
            obj.resLaser.on
            pause(0.1);
            
            
            % scan the laserpower
            if i == 1
                status.String = 'Scanning power'; drawnow;
                scanLaserPower(obj, status, Nsteps, scanPT)
            end
            
            % tune laser power to desired value
            status.String = 'Power tuning'; drawnow;
            tuneLaserPower(obj, status)
         
            
            % sweep through angles
            for k = 1:Nangles
                assert(~obj.abort_request,'User aborted.');
                theta  = rangePOL(k);
                status.String = sprintf( 'Navigating to %g (%i/%i)', theta, k, Nangles); drawnow;
                obj.rot.move(theta);
                pause4Move(obj.rot, obj.motor_move_time);
                status.String = sprintf( 'Measuring at %g (%i/%i)', theta, k, Nangles); drawnow;
               
                scanP.RData = [];
                % measure different spots
                for m = 1:length(x)
                    managers.Stages.move([x(m),y(m),zpos]);
                    managers.Stages.waitUntilStopped;
                    
                    switch obj.Detection_Type
                        case 'Spectrometer'
                            error_flag = 5;
                            while error_flag
                                assert(~obj.abort_request,'User aborted.');
                                try
                                    obj.takeSpec.run(status,managers,ax(2)) % Run spectrum
                                    spectrum = obj.takeSpec.GetData(managers.Stages,managers.Imaging);  % Get spectrum data
                                    break
                                catch err
                                    error_flag = error_flag-1;
                                    pause(10);
                                end
                            end
                            if ~error_flag
                                continue
                            end
                            obj.data.spec_wavelength = spectrum.wavelength; %
                            obj.data.spec_intensity(i,k,m,:) = spectrum.intensity;
                            %                         scanH.CData(i,:) = spectrum.intensity;
                            title(ax(2),sprintf('Spectra %i of %i',((i-1)*length(rangePOL)+k-1)*length(x)+m,length(rangeWL)*length(rangePOL)*length(x)));
                            drawnow
                        case 'APD'
                            APD_counts = obj.counter.singleShot(obj.APD_dwell)*obj.APD_dwell/1000;
                            if m == 1
                                scanP.RData(k) = APD_counts;
                            end
                            obj.data.APD_counts(i,k,m) = APD_counts;
                            drawnow
                    end
                    managers.Stages.move(startingPos);
                    managers.Stages.waitUntilStopped;
                end
            end
            
            switch obj.Detection_Type
                case 'Spectrometer'
                    scanH.CData(i,:) = squeeze(mean(obj.data.spec_intensity(i,:,1,:),2));
                    drawnow
                case 'APD'
                    scanH.YData(i) = mean(scanP.RData);
                    drawnow
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
            
            scanPW.YData(i) = nanmean(PWall(i,:))*1e3;
                    
            obj.data.rangePOL = rangePOL;
            obj.data.power_motor_pos(i) = obj.rotPower.Position;
            obj.data.laser_power(i) = nanmean(PWall(i,:)); % Save average power measured
            obj.data.laser_power_all(i,:) = PWall(i,:); % Save average power measured
            obj.data.laser_wavelength(i) = Sources.TunableLaser_invisible.c/obj.resLaser.getFrequency; % Get exact wavelength from wavemeter
            obj.resLaser.off
            
        end
        obj.data.Xcoord = x;
        obj.data.Ycoord = y;
        
        obj.data.pol_angle = rangePOL;
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

% Wait until motor stops moving, or timeout
function pause4Homing(rot,maxTime)
    t = tic;
    while (rot.Moving || ~rot.Homed) && (toc(t) < maxTime)
        drawnow
        if toc(t) > maxTime
            rot.abort;
            rot.disable
            pause(5);
            rot.enable
            pause(1);
            rot.home;
            tt=tic;
                while (rot.Moving || ~rot.Homed) && (toc(t) < maxTime)
                    drawnow
                    if toc(tt) > maxTime
                        error('Motor timed out while moving')
                    end
                end
        end
    end
end

% scan laser power
function scanLaserPower(obj, status, Nsteps, scanPT)
    PWavg_num = 25;
    position = linspace(0,2.5,Nsteps);
    delta_pos = position(2)-position(1);
    
    % always measure without OD1
    obj.PowerFlipMount.setState(2); 
    
    % get power
    for step = 1:Nsteps
        assert(~obj.abort_request,'User aborted.');
        status.String = sprintf('Scanning Laser Power %i/%i',step,Nsteps); drawnow;
        if position(step)
            obj.rotPower.step(delta_pos);
        end
        pause4Move(obj.rotPower, obj.motor_move_time*10);
        scanPT.XData(step) = obj.rotPower.Position;
        scanPT.YData(step) = get_mean_power(obj, PWavg_num)*1e3;
    end
    
    obj.data.lp_position = scanPT.XData; 
    obj.data.lp_pol = scanPT.YData*1e-3;
end

function f = sinfunc(x,xdata)
    f = (x(2)+x(1)*sin(pi/2.5*xdata+x(3)).^2);
end

% tune laser power
function tuneLaserPower(obj, status)

    PWavg_num = 25;
    
    lp = obj.data.lp_pol;
    position = obj.data.lp_position;
    [max_p, max_idx] = max(lp);
    [min_p, min_idx] = min(lp);
    orientation = sign(max_idx-min_idx);
    
    start_delta = abs( position(max_idx) - position(min_idx) ) / 5;

    % if position outside tune region: move into the region
    pos = obj.rotPower.Position;
    if ( pos > position(max_idx) && pos > position(min_idx) ) || ( pos < position(max_idx) && pos < position(min_idx) ) 
        % move to closest edge
        if  abs(pos - position(max_idx)) < abs(pos - position(min_idx))
            obj.rotPower.move(position(max_idx));
        else
            obj.rotPower.move(position(min_idx));
        end
        pause4Move(obj.rotPower, obj.motor_move_time*10);
    end
    
    power = get_mean_power(obj, PWavg_num);
    status.String = sprintf('Tuning power: %0.2f uW/%0.2f uW at pos %0.2f and OD1 %0.0f', power*1e3, obj.LaserPowerTarget*1e3, obj.rotPower.Position, obj.PowerFlipMount.getState ); drawnow;
    deltaPower = power - obj.LaserPowerTarget;
    delta = start_delta;
    end_loop = 0;
    stepped_ouside = 0;
    
    % tune power
    while( abs( deltaPower ) > obj.LaserPowerTolerance && ~end_loop )
        assert(~obj.abort_request,'User aborted.');
        
        pos = obj.rotPower.Position;

        % power too low
        if power < obj.LaserPowerTarget    
            obj.rotPower.step(delta*orientation); % make power stronger
            pause4Move(obj.rotPower, obj.motor_move_time*10);
        % power too high    
        else
            obj.rotPower.step(-delta*orientation); % weaken power
            pause4Move(obj.rotPower, obj.motor_move_time*10);
        end
        
        
        % check if we stepped outside
        if abs(pos - position(max_idx)) + abs(pos - position(min_idx)) > abs( position(max_idx)-position(min_idx) )
            % outside because too low power
            if ( orientation == 1 && pos > position(max_idx) ) || ( orientation == -1 && pos < position(max_idx) )
                % goto max
                obj.rotPower.move(position(max_idx));
                pause4Move(obj.rotPower, obj.motor_move_time*10);
                
                % check OD 1
                if obj.PowerFlipMount.getState == 1
                    obj.PowerFlipMount.setState(2); % Remove OD1 filter
                    delta = start_delta;
                    msg = 'OD1 removed.'
                else
                    if stepped_ouside
                        % cant get better, abort at maximum
                        end_loop = 1;
                        msg = 'Laser Power Tuning abort at maximum.'
                        
                    else
                        stepped_ouside = stepped_ouside + 1;
                    end
                end
            % outside because too high power
            else
                % goto min
                obj.rotPower.move(position(min_idx));
                pause4Move(obj.rotPower, obj.motor_move_time*10);
                % check OD 1
                if obj.PowerFlipMount.getState == 2
                    obj.PowerFlipMount.setState(1); % Add OD1 filter
                    delta = start_delta;
                    msg = 'OD1 added.'
                else
                    if stepped_ouside
                        % cant get better, abort at minimum
                        end_loop = 1;
                        msg = 'Laser Power Tuning abort at minima'
                    else
                        stepped_ouside = stepped_ouside + 1;
                    end
                end
            end
        end

        
        power = get_mean_power(obj, PWavg_num);
        status.String = sprintf('Tuning power: %0.2f uW/%0.2f uW at pos %0.2f and OD1 %0.0f', power*1e3, obj.LaserPowerTarget*1e3, obj.rotPower.Position, obj.PowerFlipMount.getState ); drawnow;
        old_deltaPower = deltaPower;
        deltaPower = power - obj.LaserPowerTarget;
        if sign(old_deltaPower*deltaPower) == -1
            delta = delta*0.5; % make stepsize smaller
        end
    end
    
end

% % tune laser power
% function tuneLaserPower(obj,startpos,i)
%     obj.rotPower.move(startpos);
%     pause4Move(obj.rotPower, obj.motor_move_time*10);
%     PWavg_num = 25;
%     Nsteps = 11;
%     position = linspace(0,2.5,Nsteps);
%     delta_pos = position(2)-position(1);
%     lp = zeros(1,Nsteps);
%     pos_real = zeros(1,Nsteps);
%     RepeatMeasurement = 1;
%     while (RepeatMeasurement)
%         msg='recording power'
%         for step = 1:Nsteps
%             assert(~obj.abort_request,'User aborted.');
%             if position(step)
%                 obj.rotPower.step(delta_pos);
%             end
%             pause4Move(obj.rotPower, obj.motor_move_time);
%             pos_real(step) = obj.rotPower.Position;
%             lp(step) = get_mean_power(obj, PWavg_num);
%         end
%         
%         obj.data.lp_pol(i,:)=lp;
%         obj.data.pos_pol(i,:)=pos_real;
%         
%         obj.rotPower.step(-2.5);
%         pause4Move(obj.rotPower, obj.motor_move_time*10);
%         [max_p, max_idx] = max(lp);
%         [min_p, min_idx] = min(lp);
%         
%         if max_p<obj.LaserPowerTarget && obj.PowerFlipMount.getState == 1
%             obj.PowerFlipMount.setState(2); % Remove OD1 filter
%         elseif min_p>=obj.LaserPowerTarget && obj.PowerFlipMount.getState == 2
%             obj.PowerFlipMount.setState(1); % Add OD1 filter
%         else 
%             RepeatMeasurement = 0;
%         end
%     end
%         
% %     if max_p<=obj.LaserPowerTarget && obj.PowerFlipMount.getState == 2
%     if max_p<=obj.LaserPowerTarget
%         pos = position(max_idx);
% %     elseif min_p>=obj.LaserPowerTarget && obj.PowerFlipMount.getState == 1
%     elseif min_p>=obj.LaserPowerTarget
%         pos = position(min_idx);
%     else
%         sinfunc=@(x,xdata)(x(2)+x(1)*sin(pi/2.5*xdata+x(3)).^2);
%         x0=[max(lp)-min(lp); min(lp); 1];
%         [x,resnorm] = lsqcurvefit(sinfunc,x0,position,lp);
%         x
%         deltasin=@(pos)(sinfunc(x,pos)-obj.LaserPowerTarget);
% %         pos_temp = fzero(deltasin,Nsteps*deltaAngle*0.5);
%         posrange = linspace(0,2.5,1000);
%         [pos_temp, pos_temp_idx] = min(abs(deltasin(posrange)));
%         pos = posrange(pos_temp_idx);
%     end
%     if pos
%         obj.rotPower.step(pos);
%     end
%     pause4Move(obj.rotPower, obj.motor_move_time*10);
% end

% tune laser power
% function tuneLaserPower(obj)
%     PWavg_num = 25;
%     deltaAngle = 0.125;
%     
%     initPower = get_mean_power(obj, PWavg_num);
%     initPos = obj.rotPower.Position;
%     laserpower = initPower;
%     steps = 1;
%     
%     while ( abs(laserpower-obj.LaserPowerTarget) > obj.LaserPowerTolerance)
%         deltaPower = laserpower-obj.LaserPowerTarget;
%         deltaAngle_old = deltaAngle;
%         if deltaPower < 0 % sign of power difference determines move direction
%             deltaAngle = abs(deltaAngle) * -1;
%         end
%         vz = deltaAngle*deltaAngle_old;
%         if vz < 0 % if direction changes the stepsize is halfed
%             deltaAngle = deltaAngle * 0.5;
%         end
%         
%         obj.rotPower.step(deltaAngle);
%         pause4Move(obj.rotPower, obj.motor_move_time);
%         
%         
%         laserpower = get_mean_power(obj, PWavg_num)
%         obj.rotPower.Position
%         
%         steps = steps + 1;
%         
%         if steps > 10
%             error('Laser power tuning failed')
%             obj.rotPower.move(initPos);
%             pause4Move(obj.rotPower, obj.motor_move_time);
%             break
%         end
%         
%     end
% 
% end


function NM = get_mean_power(obj,N)
    laserpower = zeros(1,N);
    for i = 1:N
        try
            laserpower(i) = obj.PM.get_power('MW');
        catch err
            warning(err.message)
            laserpower(i) = NaN;
            continue
        end
    end
    NM = nanmean(laserpower);
    return
end