function run(obj,statusH,managers,ax)
    
    obj.abort_request = false;
    status.String = 'Experiment started';
    drawnow;

    if obj.useROI
        obj.xmin = managers.Imaging.ROI(1,1);
        obj.xmax = managers.Imaging.ROI(1,2);
        obj.ymin = managers.Imaging.ROI(2,1);
        obj.ymax = managers.Imaging.ROI(2,2);
        obj.xpoints = managers.Imaging.active_module.resolution(1);
        obj.ypoints = managers.Imaging.active_module.resolution(2);
    else
        managers.Imaging.setROI([obj.xmin obj.xmax; obj.ymin obj.ymax]);
        managers.Imaging.active_module.resolution = [obj.xpoints obj.ypoints];
    end

    if strcmp(obj.scan_type,'Wavelength')
        values = obj.wavelengths;
    else
        values = obj.percents;
    end
    
    try
        for i=1:length(values)
            obj.resLaser.off;
            if strcmp(obj.scan_type,'Wavelength')
                try
                    obj.resLaser.TuneSetpoint(obj.c/values(i));
                catch
                    pause(60); % in case Msquared is doing the pump beam alignment
                    obj.resLaser.TuneSetpoint(obj.c/values(i));
                end
            else
                obj.resLaser.TunePercent(values(i));
            end

            resLaserWL = obj.resLaser.getWavelength;
            obj.resLaser.on;
            managers.Imaging.snap(true);
            if ~managers.Imaging.last_sandboxed_fn_eval_success
                error('Imaging failed');
            end

            obj.data.scans(i).image = managers.Imaging.current_image.info;
            obj.data.scans(i).Setpoint = values(i);
            obj.data.scans(i).resLaser = resLaserWL;
            obj.resLaser.off;
            
            assert(~obj.abort_request,'User aborted');
        end
    catch run_err
        if ~isempty(run_err)
            obj.resLaser.off;
            rethrow(run_err)
        end
    end

end