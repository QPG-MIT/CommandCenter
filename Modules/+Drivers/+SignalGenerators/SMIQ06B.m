classdef SMIQ06B < Drivers.SignalGenerators.SignalGenerator 
    % Matlab Object Class implementing control for SMIQ06B Signal Generator
    %
    %
    % Primary purpose of this is to control the SG
    %
    % should be named serial in an class that calls this instance
    
    % One instance controls 1 physical device. Multiple instances can exist
    % simultaneously for different hardware channels. If you need two
    % devices to work together, a new class should be designed.
    %
    % State information of the machine is stored on the SG. Can be obtained
    % using the get methods.
    
    properties
        prefs = {'comObjectInfo'};
        comObjectInfo = struct('comType','','comAddress','','comProperties','') 
        comObject;     % USB-Serial/GPIB/Prologix
    end
  
    methods(Static)
        
        function obj = instance(name)
            mlock;
            persistent Objects
            if isempty(Objects)
                Objects = Drivers.SignalGenerators.SMIQ06B.empty(1,0);
            end
            for i = 1:length(Objects)
                if isvalid(Objects(i)) && isequal(name,Objects(i).singleton_id)
                    obj = Objects(i);
                    return
                end
            end
            obj = Drivers.SignalGenerators.SMIQ06B();
            obj.singleton_id = name;
            Objects(end+1) = obj;
        end
        
    end
    
    methods(Access=private)
        function [obj] = SMIQ06B()
            obj.loadPrefs;
            display('setting comInfo for SMIQ06B.')
            if isempty(obj.comObjectInfo.comType)&& isempty(obj.comObjectInfo.comAddress)&& isempty(obj.comObjectInfo.comProperties)
                 %first time connecting should run the helperfunction
                %Connect_Device to establish your connection
                [obj.comObject,obj.comObjectInfo.comType,obj.comObjectInfo.comAddress,obj.comObjectInfo.comProperties] = Connect_Device;
            else
                try
                    %this is used for connecting every time after the first
                    %time
                    [obj.comObject,obj.comObjectInfo.comType,obj.comObjectInfo.comAddress,obj.comObjectInfo.comProperties] = ...
                        Connect_Device(obj.comObjectInfo.comType,obj.comObjectInfo.comAddress,obj.comObjectInfo.comProperties);
                catch
                    %this is only called if you change a device property
                    %after the intiial connection (ex: change GPIB
                    %address). This allows you to establish a new
                    %connection.
                    [obj.comObject,obj.comObjectInfo.comType,obj.comObjectInfo.comAddress,obj.comObjectInfo.comProperties] ...
                        = Connect_Device;
                end
            end
            fopen(obj.comObject);
            obj.reset;
        end
    end
    
    methods(Access=private)
        function DeleteListFreq(obj)
            string  = sprintf('LIST:DELete:FREQ');
            obj.writeOnly(string);
        end
        
        function DeleteListPower(obj)
            string  = sprintf('LIST:DELete:POWer');
            obj.writeOnly(string);
        end
        
        function  ListLearn(obj)
            string = sprintf('SOURCE:LIST:LEARN');
            obj.writeOnly(string);
        end
        
        function writeOnly(obj,string)
            fprintf(obj.comObject,string);
        end
        
        function [output] = writeRead(obj,string)
            output = query(obj.comObject,string);
        end
        
        function  setListTrig(obj,ListTrig)
            string = sprintf('TRIGGER:LIST:SOURCE %s',ListTrig);
            obj.writeOnly(string);
        end
        
        function  setFreqMode(obj,FreqMode)
            switch FreqMode
                case {'FIX','CW'},
                    obj.writeOnly('SOUR:FREQ:MODE FIX');
                case {'LIST'}
                    obj.writeOnly('SOURCE:FREQUENCY:MODE LIST');
                    obj.writeOnly('SOURCE:LIST:MODE STEP');     %Only STEP modes can be used
                otherwise
                    warning('No frequency mode was set');
            end
        end
        
        function  setPowerMode(obj,PowerMode)
            string = sprintf('SOURce:POWer:MODE %s',PowerMode);
            obj.writeOnly(string);
        end
    end
    methods
        
        
        function setUnitPower(obj)
            obj.writeOnly('UNIT:POWER DBM');
        end
        
        function  setFreqCW(obj,Freq)
            string = sprintf(':FREQuency:FIXed %f',Freq);
            obj.writeOnly(string);
        end
        
        function  setPowerCW(obj,Power)
            string = sprintf(':POWER:LEVEL:IMMEDIATE:AMPLITUDE %f DBM',Power);
            obj.writeOnly(string);
        end
        
        function  setFreqList(obj,FreqList)
            
            obj.DeleteListFreq();
            
            obj.writeOnly('*WAI');
            
            NumberOfPoints = length(FreqList);
            
            clear string;
            
            if NumberOfPoints > 0
                for i = 1:NumberOfPoints
                    if i == 1
                        string = sprintf('%f',FreqList(i));
                    else
                        string = [string sprintf(',%f',FreqList(i))];
                    end
                end
                string = [':LIST:FREQUENCY ',string];
                obj.writeOnly(string);
            end
            obj.writeOnly('*WAI');
        end
        
        function  setPowerList(obj,PowerList)
            
            obj.DeleteListPower();
            
            obj.writeOnly('*WAI');
            NumberOfPoints = length(PowerList);
            
            clear string;
            
            if NumberOfPoints > 0
                for i = 1:NumberOfPoints
                    if i == 1
                        string = sprintf('%f',PowerList(i));
                    else
                        string = [string sprintf(',%f',PowerList(i))];
                    end
                end
                string = ['SOURCE:LIST:POWER ',string];
                obj.writeOnly(string);
            end
            obj.writeOnly('*WAI');
        end
        
        %%
        
        function UnitPower = getUnitPower(obj)
            UnitPower = obj.writeRead('UNIT:POWER?');
            UnitPower = strrep(UnitPower,sprintf('\n'),''); %remove excess carriage returns
        end
        
        function  [Freq]=getFreqCW(obj)
            string = sprintf('FREQ?');
            s = obj.writeRead(string);
            Freq = str2num(s);
        end
        
        function  [Power]=getPowerCW(obj)
            string = sprintf('POW?');
            s = obj.writeRead(string);
            Power = str2num(s);
        end
        
        function  [FreqMode]=getFreqMode(obj)
            string = sprintf('FREQ:MODE?');
            s = obj.writeRead(string);
            FreqMode = strrep(s,sprintf('\n'),''); %remove excess carriage returns;
        end
        
        function  [PowerMode]=getPowerMode(obj)
            string = sprintf('POWer:MODE?');
            s = obj.writeRead(string);
            PowerMode = strrep(s,sprintf('\n'),''); %remove excess carriage returns;
        end
        
        function  [FreqList]=getFreqList(obj)
            string = sprintf('LIST:FREQ?');
            s = obj.writeRead(string);
            FreqList = str2num(s);
        end
        
        function  [PowerList]=getPowerList(obj)
            string = sprintf('LIST:POWer?');
            s = obj.writeRead(string);
            PowerList = str2num(s);
        end
        
        function  [MWstate]=getMWstate(obj)
            string = sprintf('OUTPUT:STATE?');
            s = obj.writeRead(string);
            if strcmp(strrep(s,sprintf('\n'),''),'1')
                MWstate = 'On';
            else
                MWstate = 'Off';
            end
        end
        
        function program_list(obj,freq_list,power_list)
            obj.reset;
            obj.setFreqMode('CW');
            obj.setPowerMode('CW');
            obj.select_list('LIST1')
            obj.setFreqList(freq_list);
            obj.setPowerList(power_list);
            obj.ListLearn;
            obj.setFreqMode('LIST');
            obj.setPowerMode('LIST');
            obj.setListTrig('EXT');
            obj.off;
        end
        
        function select_list(obj,listname)
            assert(ischar(listname),'SMIQ list name must be a string.')
            obj.writeOnly(['SOUR:LIST:SEL ''' listname ''''])
        end
        
        %% 
        function  on(obj)
            string = sprintf(':OUTPUT:STATE 1');
            obj.writeOnly(string);
        end
        
        function  off(obj)
            string = sprintf(':OUTPUT:STATE 0');
            obj.writeOnly(string);
        end
        
        function delete(obj)
            obj.reset;
            fclose(obj.comObject);
            delete(obj.comObject);
        end
        
        function  reset(obj)
            string = sprintf('*RST');
            obj.writeOnly(string);
        end
        
    end
end