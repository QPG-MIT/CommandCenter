function run( obj,status,managers,ax )
    % Main run method (callback for CC run button)
    obj.abort_request = false;
    status.String = 'Experiment started';
    drawnow;
    % Edit here down (save data to obj.data)
    % Tips:
    % - If using a loop, it is good practice to call:
    %     drawnow; assert(~obj.abort_request,'User aborted.');
    %     as frequently as possible
    % - try/catch/end statements useful for cleaning up
    % - You can get a figure-like object (to create subplots) by:
    %     panel = ax.Parent; delete(ax);
    %     ax(1) = subplot(1,2,1,'parent',panel);
    % - drawnow can be used to update status box message and any plots

    % Edit this to include meta data for this experimental run (saved in obj.GetData)
    obj.meta.prefs = obj.prefs2struct;
    obj.meta.position = managers.Stages.position; % Save current stage position (x,y,z);
    obj.meta.angles = obj.angles; %Angles corresponding to each spectrum

    % Check that rot is not empty and valid
    assert(~isempty(obj.rot) && isvalid(obj.rot),'Motor SN must be a valid number.')

    try
        % Instantiate driver for the rotation mount
        rot.home()
        while rot.isMoving()
            drawnow
        end

        % Sweep through polarisation and get spectra
        for theta = obj.angle_list
            obj.rot.move(theta)
            while rot.isMoving()
                drawnow
            end

            RunExperiment(obj, managers, obj.spec_experiment, theta, ax)
            %obj.data.angle(theta) = obj.spec_experiment.GetData
            tempDat = obj.spec_experiment.GetData;
            obj.data.angle(theta).wavelength = tempDat.wavelength;
            obj.data.angle(theta).intensity = tempDat.intensity;
            drawnow; assert(~obj.abort_request,'User aborted');
        end

        obj.meta.spec_meta = obj.spec_experiment.meta; %Get meta data from spectrum experiment



    catch err
    end
    % CLEAN UP CODE %
    if exist('err','var')
        % HANDLE ERROR CODE %
        rethrow(err)
    end
end

function RunExperiment(obj,managers,experiment,site_index,ax)
    [abortBox,abortH] = ExperimentManager.abortBox(class(experiment),@(~,~)obj.abort);
    try
        drawnow; assert(~obj.abort_request,'User aborted');
        if ~isempty(experiment.path) %if path defined, select path
            managers.Path.select_path(experiment.path);
        end
        obj.current_experiment = experiment;
        cla(ax,'reset');
        experiment.run(abortBox,managers,ax);
        obj.current_experiment = [];
    catch exp_err
        obj.data.sites(site_index).experiments(end).err = exp_err;
        delete(abortH);
        rethrow(exp_err)
    end
    delete(abortH);
end