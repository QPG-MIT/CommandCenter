function view(obj)
if ~isempty(obj.GUI)
    figure(obj.GUI.fig)
    return
end
GUI.fig = figure('IntegerHandle','off',...
    'name',sprintf('NIDAQ: %s',obj.DeviceChannel),...
    'visible','off','HandleVisibility','Callback',...
    'units','inches',...
    'CloseRequestFcn',@obj.close,...
    'menu','none');
pos = GUI.fig.Position;
pos(3:4) = [4 2];
GUI.fig.Position = pos;

%% Setup item lists
% Tasks
p = uipanel(GUI.fig,'title','Tasks',...
    'Position',[0 0 1/3 1]);
GUI.tasks = uicontrol(p,'Style','listbox',...
    'units','normalized',...
    'tooltipstring',sprintf('Started - green\nStopped - black\nAborted - red\nOtherwise - name: status'),...
    'Position',[0.01 0.01 1 0.97],...
    'enable','inactive',...
    'tooltipstring','Name',...
    'value',[],'max',Inf,'min',0);
% In Lines
p = uipanel(GUI.fig,'title','In Lines',...
    'Position',[1/3 0 1/3 1]);
GUI.InLines = uicontrol(p,'Style','listbox',...
    'units','normalized',...
    'tooltipstring','Name: Hardware Line',...
    'Position',[0.01 0.01+0.12 1 0.97-0.12]);
uicontrol(p,'Style','Pushbutton',...
    'units','normalized',...
    'Position',[0.05 0.025 0.45 0.1],...
    'String','New',...
    'Callback',@obj.addInLine_Callback)
uicontrol(p,'Style','Pushbutton',...
    'units','normalized',...
    'Position',[0.5 0.025 0.45 0.1],...
    'String','Remove',...
    'Callback',@obj.removeInLine_Callback)
% Out Lines
p = uipanel(GUI.fig,'title','Out Lines',...
    'Position',[2/3 0 1/3 1]);
GUI.OutLines = uicontrol(p,'Style','listbox',...
    'tooltipstring',sprintf('Name: Hardware Line (State)\nState is in volts for analog.'),...
    'units','normalized',...
    'Position',[0.01 0.01+0.12 1 0.97-0.12]);
uicontrol(p,'Style','Pushbutton',...
    'units','normalized',...
    'Position',[0.05 0.025 0.45 0.1],...
    'String','New',...
    'Callback',@obj.addOutLine_Callback)
uicontrol(p,'Style','Pushbutton',...
    'units','normalized',...
    'Position',[0.5 0.025 0.45 0.1],...
    'String','Remove',...
    'Callback',@obj.removeOutLine_Callback)

obj.GUI = GUI;
obj.GUI.listeners(1)=addlistener(obj,'InLines','PostSet',@obj.RefreshLines);
obj.GUI.listeners(1)=addlistener(obj,'OutLines','PostSet',@obj.RefreshLines);
obj.GUI.listeners(2)=addlistener(obj,'Tasks','PostSet',@obj.RefreshTasks);
GUI.fig.Visible = 'on';
obj.RefreshTasks;
obj.RefreshLines;
for i = 1:numel(obj.OutLines)
    line = obj.OutLines(i);
    line.niListener = addlistener(line,'state','PostSet',@obj.UpdateLine);
end
for i = 1:numel(obj.Tasks)
    task = obj.Tasks(i);
    task.niListener = addlistener(task,'status','PostSet',@obj.UpdateTask);
end

drawnow
