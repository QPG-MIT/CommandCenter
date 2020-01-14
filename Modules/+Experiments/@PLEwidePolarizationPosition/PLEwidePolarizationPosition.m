classdef PLEwidePolarizationPosition < Modules.Experiment
    % PLEwidePolarization Description of experiment
    % Useful to list any dependencies here too
    
    properties(SetObservable,AbortSet)
        
        resLaser = Modules.Source.empty(1,0); % Call EMM or solstis
        takeSpec = Experiments.Spectrum.instance; % Call Spectrometer
        ScanRange = 'linspace(710,720,11)'; % %eval(ScanRange) will define Scanning range (Range)
        Detection_Type = {'APD','Spectrometer'}; 
        APD_line = 'APD1'; % APD line name (for counter)
        APD_dwell = 100; % APD dwell time [ms]
%         Counter = Modules.Driver.Counter.instance; % Call Counter for APD
%         hwserver_ip = Experiments.PLE_WideScan_SPEC.no_server;
        PolarizationAngleRange = 'linspace(0,180,11)'; % 
        LaserPowerTarget = 1; % laser power targe in [uW]
        LaserPowerTolerance = 0.01; % laser power tolerance in [uW]
        Xcoords = '[0]';
        Ycoords = '[0]';
%         motor_SN = Drivers.APTMotor.getAvailMotors; % Serial number for the rotation mount, to be used to create a driver for the rotation mount. Must be connected through APT Config first.
%         power_motor_SN = Drivers.APTMotor.getAvailMotors; % Serial number for the rotation mount, to be used to create a driver for the rotation mount. Must be connected through APT Config first.
        motor_SN = '27251915'; % Serial number for the rotation mount, to be used to create a driver for the rotation mount. Must be connected through APT Config first.
        power_motor_SN = '83847313';
%         power_flipmount_SN = '37000086';
        motor_move_time = 30;  % Maximum time allowed for motor to move between positions
        motor_home_time = 120; % Maximum time allowed for the motor to home itself
    end
    properties
        prefs = {'resLaser','takeSpec','ScanRange','PolarizationAngleRange','Xcoords','Ycoords','LaserPowerTarget','LaserPowerTolerance','Detection_Type','APD_line','APD_dwell','motor_SN','power_motor_SN'};  % String representation of desired prefs
        %show_prefs = {'resLaser','takeSpec','ScanRange','PolarizationAngleRange','Detection_Type','APD_line','APD_dwell'};   % Use for ordering and/or selecting which prefs to show in GUI
        %readonly_prefs = {}; % CC will leave these as disabled in GUI (if in prefs/show_prefs)
        PM;
        counter;                % APD counter
        rot;
        rotPower;
        PowerFlipMount;
    end
    properties(SetAccess=private,Hidden)
        % Internal properties that should not be accessible by command line
        % Advanced users should feel free to alter these properties (keep in mind methods: abort, GetData)
        data = [] % Useful for saving data from run method
        meta = [] % Useful to store meta data in run method
        abort_request = false; % Flag that will be set to true upon abort. Use in run method!
    end
    properties(Constant,Hidden)
        no_server = 'No Server';  % Message when not connected
    end

    methods(Static)
        % Static instance method is how to call this experiment
        % This is a separate file
        obj = instance()
    end
    methods(Access=private)
        function obj = PLEwidePolarizationPosition()
            % Constructor (should not be accessible to command line!)
            obj.loadPrefs; % Load prefs specified as obj.prefs
%             obj.ni = Drivers.NIDAQ.dev.instance('dev1');
            obj.PM = Drivers.PM100.instance;
            obj.counter = Drivers.Counter.instance(obj.APD_line,'CounterSync');
            obj.PowerFlipMount = Drivers.APTFilterFlipper.instance(37000086);
            obj.rot = Drivers.APTMotor.instance(27251915, [0 360]);
            obj.rotPower = Drivers.APTMotor.instance(83847313, [-180 180]);
        end
    end

    methods
        run(obj,status,managers,ax) % Main run method in separate file

        function abort(obj)
            % Callback for when user presses abort in CC
            obj.abort_request = true;
        end

        function dat = GetData(obj,stageManager,imagingManager) 
            % Callback for saving methods (in CommandCenter)
%             meta = stageManager.position;
            dat.data = obj.data;
            dat.meta = obj.meta;
        end

        % Set methods allow validating property/pref set values
        function set.resLaser(obj,val)
            if isempty(val)
                obj.resLaser = val;
                return
            end
            h = superclasses(val);
            assert(ismember('Sources.TunableLaser_invisible',h),'Laser must be tunable')
            obj.resLaser = val;
        end
        function set.ScanRange(obj,val)
            t = eval(val);
            assert(~isempty(t),'ScanRange is empty')
            assert(isnumeric(t),'Value must be numeric')
            obj.ScanRange = val;
        end
        function set.PolarizationAngleRange(obj,val)
            t = eval(val);
            assert(~isempty(t),'ScanRange is empty')
            assert(isnumeric(t),'Value must be numeric')
            obj.PolarizationAngleRange = val;
        end
        function set.Xcoords(obj,val)
            t = eval(val);
            assert(~isempty(t),'Xcoords is empty')
            assert(isnumeric(t),'Value must be numeric')
            obj.Xcoords = val;
        end
        function set.Ycoords(obj,val)
            t = eval(val);
            assert(~isempty(t),'Ycoords is empty')
            assert(isnumeric(t),'Value must be numeric')
            obj.Ycoords = val;
        end
%         
%         function set.motor_SN(obj,val)
% %             val = '27251915';
%             if ~isempty(val)
%                 val_as_double = str2double(val); % must be double to instantiate motor
%                 assert(~isnan(val_as_double),'Motor SN must be a valid number.')
% 
%                 Handle proper deleting of smotor driver object
%                 try
%                     delete(obj.rot); % Either motor obj or empty
%                     obj.rot = [];
%                 catch
%                 end
%                 
%                 obj.motor_SN = val;
%                 if val_as_double == 0
%                     %Leave obj.rot empty if no serial number selected
%                     return % Short circuit
%                 end
% 
%                 % Add new motor
%                 try
%                     obj.rot = Drivers.APTMotor.instance(val_as_double, [0 360]);
%                 catch
%                     obj.rot = [];
%                 end 
%             end
%         end
%         
%         function set.power_motor_SN(obj,val)
% %             val = '83847313';
%             if ~isempty(val)
%                 val_as_double = str2double(val); % must be double to instantiate motor
%                 assert(~isnan(val_as_double),'Motor SN must be a valid number.')
% 
%                 % Handle proper deleting of smotor driver object
%                 try
%                     delete(obj.rotPower); % Either motor obj or empty
%                     obj.rotPower = [];
%                 catch 
%                 end
%                 obj.power_motor_SN = val;
%                 if val_as_double == 0
%                     %Leave obj.rotPower empty if no serial number selected
%                     return % Short circuit
%                 end
% 
%                 % Add new motor
%                 try
%                     obj.rotPower = Drivers.APTMotor.instance(val_as_double, [0 360]);
%                 catch
%                     obj.rotPower = [];
%                 end
%         
%             end
%         end
%         function set.power_flipmount_SN(obj,val)
% %             val = '83847313';
%             if ~isempty(val)
%                 val_as_double = str2double(val) % must be double to instantiate motor
%                 assert(~isnan(val_as_double),'Flipmount SN must be a valid number.')
% 
%                 % Handle proper deleting of smotor driver object
%                 obj.power_flipmount_SN = val;
%                 if val_as_double == 0
%                     %Leave obj.rotPower empty if no serial number selected
%                     return % Short circuit
%                 end
% 
%                 % Add new flipper
%                 try
%                     obj.PowerFlipMount = Drivers.APTFilterFlipper.instance(val_as_double);
%                 catch
%                     obj.PowerFlipMount = [];
%                 end
%         
%             end
%         end
        
        function delete(obj)
            delete(obj.rot);
            delete(obj.rotPower);
        end
    end
end
