function varargout = wfanalyze(data,varargin)
%ANALYZE Examine data from all sites
%   Inputs:
%       data: data produced by GetData method
%       [Analysis]: An analysis struct produced by this function
%       [FitType]: "gauss", "lorentz", or "voigt" (default "gauss")
%       [inds]: array of indices to mask full dataset (default: all data).
%          This will also filter analysis if provided.
%       [viewonly]: do not begin uifitpeaks on the axes
%       [new]: arrow navigation will go to nearest new site (e.g. continued = 0)
%       [preanalyze]: if true, will fit all data before opening UI
%       [block]: (false) Calls uiwait internally and will return analysis
%          when user closes figure.
%   Outputs:
%       (fig): Figure handle
%       (analysis): analysis struct (use with block=true).
%   Interactivity:
%       click on a spot to see corresponding data
%       [alt+] left/right arrows to change site site_index. The alt is only
%           necessary if viewonly=false, which is default.
%       Right clicking on axis will allow you to choose lorentz/gauss for that axis
%   Tag info for data:
%       Axes from left -> right:
%           'SpatialImageAx', 'SpectraAx', 'OpenLoopAx', 'ClosedLoopAx'
%       Image in 'SpatialImageAx': 'SpatialImage'; Scatter plot: 'sites'
%       All lines in 'SpectraAx': 'Spectra'
%       All errorfill children (line & patch) in 'OpenLoopAx': 'OpenLoop'
%       All errorfill children (line & patch) in 'ClosedLoopAx': 'ClosedLoop'
%   When viewonly = false, main keyboard functionality goes to UIFITPEAKS
%       Click on circle node to select it, use arrows to change its
%       location (and corresponding guess value).
%       ctl+arrows allow fine control.
%       [Shift+] Tab changes selected point
%   Analysis data:
%     Nx3 struct array with fields: (N is number of sites, 3 corresponds to experiments)
%       amplitudes - Nx1 double
%       widths - Nx1 double (all FWHM)
%       locations - Nx1 double
%       etas - Nx1 double
%       background - 1x1 double
%       fit - cfit object or empty if no peaks found
%       index - index into data.sites. If NaN, this wasn't analyzed
%   Can export/save from file menu
%   Will be prompted to export analysis data to base workspace upon closing
%   figure if no export since last analysis data update_all (NOTE this is only
%   saved when switching sites).
%       This will not overwrite previously exported data sets unless
%       specified.

p = inputParser();
addParameter(p,'Analysis',[],@isstruct);
addParameter(p,'FitType','gauss',@(x)any(validatestring(x,{'gauss','lorentz','voigt'})));
addParameter(p,'inds',1:length(data.data.sites),@(n)validateattributes(n,{'numeric'},{'vector'}));
addParameter(p,'viewonly',false,@islogical);
addParameter(p,'preanalyze',false,@islogical);
addParameter(p,'new',false,@islogical);
addParameter(p,'block',false,@islogical);
parse(p,varargin{:});
if p.Results.preanalyze && p.Results.viewonly
    error('SpecSlowScan analysis cannot run preanalysis in viewonly mode.')
end

prefs = data.meta.prefs; %load prefs from exp
FullData = data;
data = FullData.data; %reshape hack
im = data.image.image; %initialize image
sites = data.sites(p.Results.inds); %option to pick out which sites to look at


fig = figure('name',mfilename,'numbertitle','off'); %make a figure
fig.Position(3) = fig.Position(3)*2; %make it biggg
file_menu = findall(fig,'tag','figMenuFile'); %function calls to file menus
uimenu(file_menu,'Text','Go to Index','callback',@go_to,'separator','on');
uimenu(file_menu,'Text','Save Analysis','callback',@save_data);
uimenu(file_menu,'Text','Export Analysis','callback',@export_data);
uimenu(file_menu,'Text','Diagnostic Plot','callback',@open_diagnostic);
bg(1) = uipanel(fig,'units','normalized','position',[0   0 1/2 1],'BorderType','none');
bg(2) = uipanel(fig,'units','normalized','position',[1/2 0 1/2 1],'BorderType','none');
splitPan(1) = Base.SplitPanel(bg(1),bg(2),'horizontal'); 
set(splitPan(1).dividerH,'BorderType','etchedin') %draw a panel on the left side, where bg(2) is

% Add some help tooltips
selector(2).Tooltip = 'Dashed line corresponds to setpoint of scan at 50%.';

ax = axes('parent',bg(1),'tag','SpatialImageAx'); %make leftmost axes, where im goes
hold(ax,'on');

%plot image 
if ~isempty(im)
    imagesc(ax,im.ROI(1,:),im.ROI(2,:),im.image,'tag','SpatialImage');
end

positions = reshape([sites.position],length(data.sites(1).position),[]);
sc = scatter(positions(1,:),positions(2,:),'ButtonDownFcn',@selectSite,...
    'MarkerEdgeAlpha',0.3,'tag','sites');
sc.UserData.fig = fig;
pos = scatter(NaN,NaN,'r');
xlabel(ax,'X Position (um)');
ylabel(ax,'Y Position (um)');
colormap(fig,'gray');
axis(ax,'image');
set(ax,'ydir','normal');
hold(ax,'off');


ax(2) = axes('parent',bg(2),'tag','WfPleAx'); hold(ax(2),'on');
addlistener(ax(2),'XLim','PostSet',@xlim_changed);

% Constants and large structures go here
block = p.Results.block;
n = length(sites);
filter_new = p.Results.new;
viewonly = p.Results.viewonly;
preanalyze = p.Results.preanalyze;
FitType = p.Results.FitType;
wavenm_range = 299792./prefs.freq_range; % Used when plotting
inds = p.Results.inds;
AmplitudeSensitivity = 1;
update_exp = {@update_wfPle}; % Easily index by exp_id
colors = lines(7);

% Frequently updated and small stuff here
site_index = 1;
busy = false;
new_data = false;
if isstruct(p.Results.Analysis)
    analysis = p.Results.Analysis;
    % Backwards compatibility
    if ~isfield(analysis,'sites')
        analysis = struct('sites',analysis);
        analysis.nm2THz = [];
        analysis.gof = [];
        warning('Found old format of analysis; updated to new format.')
        new_data = true;
    end
    if ~isfield(analysis.sites,'redo')
        for isite = 1:size(analysis,1)
            for jexp = 1:size(analysis,2)
                analysis.sites(isite,jexp).redo = false;
            end
        end
        new_data = true;
        warning('Added redo flag to loaded analysis.')
    end
    if ~isfield(analysis.sites,'ignore')
        for isite = 1:size(analysis,1)
            for jexp = 1:size(analysis,2)
                analysis.sites(isite,jexp).ignore = [];
            end
        end
        new_data = true;
        warning('Added ignore flag to loaded analysis.')
    end
    if size(analysis.sites,2) == 3
        for ii = 1:size(sites,1)
            analysis.sites(ii,4).amplitudes = NaN;
            analysis.sites(ii,4).widths = NaN;
            analysis.sites(ii,4).locations = NaN;
            analysis.sites(ii,4).etas = NaN;
            analysis.sites(ii,4).background = NaN;
            analysis.sites(ii,4).index = NaN;
            analysis.sites(ii,4).redo = false;
        end
        new_data = true;
        warning('Added 4th column to analysis.sites')
    end
    % Filter with inds
    analysis.sites = analysis.sites(p.Results.inds,:);
else
    analysis.nm2THz = [];
    analysis.gof = [];
    analysis.sites = struct(...
        'fit',cell(n,4),...
        'amplitudes',NaN,...
        'widths',NaN,...
        'locations',NaN,...
        'etas',NaN,...
        'background',NaN,...
        'index',NaN,... % Index into sites
        'redo',false,...
        'ignore',[]);     % indices of experiments in the fit
end
% Can now link closereq function (prevents losing unsaved data)
fig.CloseRequestFcn = @closereq;

% Evaluate analysis for user by cycling through sites
if preanalyze
    err = containers.Map;
    site_index = 1;
    progbar = waitbar(site_index/n,'','Name','Analyzing all data','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(progbar,'canceling',0);
    try % first site has to happen out of loop, needs its own try-catch
        update_all() % call directly to bypass save_state() the first time around
    catch curr_err % save site and error message
        errtxt = getReport(curr_err);
        err(errtxt) = site_index; % Must be first error, so no need to check key
    end
    for curr_index=2:n
        try
            waitbar(curr_index/n,progbar,sprintf('Analyzing site %i/%i',curr_index,n));
            drawnow limitrate
            changeSite(curr_index);
            if getappdata(progbar,'canceling')
                break
            end
        catch curr_err
            errtxt = getReport(curr_err);
            if err.isKey(errtxt)
                err(errtxt) = [err(errtxt), site_index];
            else
                err(errtxt) = site_index;
            end
            continue
        end
    end
    delete(progbar)
    save_state()
    if err.Count % number of errors > 0
        allsites = []; % Build up all sites to list at top
        errtxt = err.keys(); % Cell array of all unique error reports
        for j = 1:err.Count
            errsites = err(errtxt{j});
            allsites = [allsites errsites];
            errsites = ['Sites ' sprintf('%g, ',err(errtxt{j}))]; % Format
            errsites = errsites(1:end-2); % remove extra ", "
            errtxt{j} = [errsites newline errtxt{j}]; % Format
        end
        allsites = sort(allsites);
        allsites = sprintf('%g, ', allsites); % Format
        allsites = allsites(1:end-2); % remove extra ", "
        delim = '---------------------------------------';
        msg = sprintf('Errors on sites %s:\n%s\n%s',...
            allsites, delim, strjoin(errtxt,...
            [newline delim newline]));
        warning(msg);
        errordlg(sprintf('Errors on sites %s. See warning in console.',allsites));
    end
end

% Link UI control
set([fig, selector],'KeyPressFcn',@cycleSite);
update_all(); % Bypass changeSite since we have no previous site

block = p.Results.block;
if block
    uiwait(fig);
end

if nargout
    varargout = {fig,analysis};
end

    function open_diagnostic(varargin)
        save_state();
        try
            [nm2THz,gof] = Experiments.AutoExperiment.SpecSlowScan.diagnostic(FullData,analysis.sites);
            if ~isequal(nm2THz,0)
                answer = questdlg('New winspec calibration fit using analyzed peaks in data; add to analysis?','WinSpec Calibration Fit','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    analysis.nm2THz = nm2THz;
                    analysis.gof = gof;
                    new_data = true;
                end
            end
        catch err
            errordlg(getReport(err,'extended','hyperlinks','off'));
        end
    end

    function export_data(varargin)
        if nargin < 1 || ~isa(fig,'matlab.ui.Figure')
            [~,fig] = gcbo;
        end
        save_state();
        if ~isempty(analysis)
            var_name = 'analysis';
            i = 1;
            while evalin('base', sprintf('exist(''%s'',''var'') == 1',var_name))
                i = i + 1;
                var_name = sprintf('%s%i','analysis',i);
            end
            if i > 1
                answer = questdlg(sprintf('Would you like to export "analysis" data to workspace as new variable "%s" or overwrite existing "analysis"?',...
                    var_name),'Export','Overwrite','New Variable','No','Overwrite');
                if strcmp(answer,'Overwrite')
                    answer = 'Yes';
                    var_name = 'analysis';
                end
            else
                answer = questdlg(sprintf('Would you like to export "analysis" data to workspace as new variable "%s"?',var_name),...
                    'Export','Yes','No','Yes');
            end
            if strcmp(answer,'Yes')
                assignin('base',var_name,analysis)
            end
        end
        new_data = false;
    end
    function save_data(varargin)
        save_state();
        last = '';
        namespace = Base.Module.get_namespace(mfilename('class'));
        if ispref(namespace,'last_save')
            last = getpref(namespace,'last_save');
        end
        [file,path] = uiputfile('*.mat','Save Analysis',last);
        if ~isequal(file,0)
            setpref(namespace,'last_save',path);
            [~,~,ext] = fileparts(file);
            if isempty(ext) % Add extension if not specified
                file = [file '.mat'];
            end
            save(fullfile(path,file),'-struct','analysis');
        end
        new_data = false;
    end

    function closereq(~,~)
        save_state();
        if ~block % If block; we are returning data; assume no save expected
            % Export data to workspace if analysis exists
            try
                if new_data
                    export_data(fig);
                end
            catch err
                delete(fig)
                rethrow(err);
            end
        end
        delete(fig)
    end

    function changeSite(new_index)
        % Only function allowed to update site_index
        if busy; return; end
        % Save the current analysis before moving to next site
        save_state();
        site_index = new_index;
        update_all();
    end

    function selectSite(sc,eventdata)
        if eventdata.Button == 1
            [~,D] = knnsearch(eventdata.IntersectionPoint(1:2),[sc.XData; sc.YData]','K',1);
            [~,ind] = min(D);
            changeSite(ind);
        end
    end

    function go_to(~,~)
        site = inputdlg(sprintf('Jump to site (between 1 and %i):',n),mfilename,1,{num2str(n)});
        if ~isempty(site)
            site_num = str2double(site{1});
            if ~isnan(site_num) && site_num <= n && site_num > 0
                changeSite(site_num);
            else
                errordlg(sprintf('"%s" is not a number between 1 and %i.',site{1},n),mfilename);
            end
        end
    end

    function cycleSite(~,eventdata)
        ind = site_index;
        for i = 1:n % Just go through sites once
            switch eventdata.Key
                case 'leftarrow'
                    direction = -1;
                case 'rightarrow'
                    direction = 1;
                otherwise % Ignore anything else
                    return
            end
            ind = mod(ind-1+direction,n)+1;
            if filter_new && ~any([sites(ind).experiments.continued]==0)
                continue
            end
            changeSite(ind);
            return
        end
        errordlg('No new sites found; try relaunching without new flag set to true.')
    end
%% Update UI methods
    function prepUI(ax,selector)
        set(selector,'Data',cell(0,10)); % Reset selector
        cla(ax,'reset'); hold(ax,'on');
    end
    function update_all()
        update_im();
        update_wfple();
    end
    function update_im()
        % Update spectrometer data (ax(1))
        site = sites(site_index);
        ax(1).Title.String = sprintf('Site %i/%i',site_index,n);
        set(pos,'xdata',site.position(1),'ydata',site.position(2));
    end
    
    function update_wfple()
        % Update PLE closed (analysis.sites(:,1), selector(2), ax(2))
        if busy; error('Busy!'); end
        busy = true;
        try
            site = sites(site_index);
            prepUI(ax(2),selector(2));
            exp_inds = fliplr(find(strcmp('Experiments.SlowScan.Closed',{site.experiments.name})));
            for i = 1:length(exp_inds)
                experiment = site.experiments(exp_inds(i));
                if ~isempty(experiment.data) &&  ~any(exp_inds(i) == analysis.sites(site_index,3).ignore)
                    errorfill(experiment.data.data.freqs_measured,...
                        experiment.data.data.sumCounts,...
                        experiment.data.data.stdCounts*sqrt(experiment.prefs.samples),...
                        'parent',ax(2),'tag','wfple','color',colors(mod(i-1,size(colors,1))+1,:));
                    formatSelector(selector(2),experiment,exp_inds(i),3,site_index,colors(mod(i-1,size(colors,1))+1,:));
                else
                    formatSelector(selector(2),experiment,exp_inds(i),3,site_index);
                end
            end
            ax(2).Title.String = 'Widefield PLE Scan';
            ax(2).XLabel.String = 'Frequency (THz)';
            ax(2).YLabel.String = 'Counts';
            if ~viewonly && ~isempty(findall(ax(2),'type','line'))
                attach_uifitpeaks(ax(4),analysis.sites(site_index,3),...
                    'AmplitudeSensitivity',AmplitudeSensitivity);
            end
        catch err
            busy = false;
            rethrow(err);
        end
        busy = false;
    end
    
    function formatSelector(selectorH,experiment,i,exp_ind,site_ind,rgb)
        if nargin < 6
            color = '';
        else
            hex = sprintf('#%02X%02X%02X',round(rgb*255));
            color = sprintf('<html><font color="%s">&#9724;</font></html>',hex);
        end
        date = datestr(experiment.tstart);
        if ~isempty(experiment.err)
            date = sprintf('<html><font color="red">%s</font></html>',date);
        elseif experiment.completed && ~experiment.skipped
            date = sprintf('<html><font color="green">%s</font></html>',date);
        end
        if isempty(experiment.tstop) % Errors cause empty
            duration = '-';
        else
            duration = char(experiment.tstop - experiment.tstart);
        end
        displayed = false;
        if ~any(i==analysis.sites(site_ind,exp_ind).ignore)
            displayed = true;
        end
        % Analysis redo should be flexible
        analysis_redo = false;
        if ~isempty(analysis.sites(site_ind,exp_ind).redo)
            analysis_redo = analysis.sites(site_ind,exp_ind).redo;
        end
        selectorH.Data(end+1,:) = {displayed,color, i,...
            date,...
            experiment.continued,...
            analysis_redo,...
            duration,...
            experiment.skipped,...
            experiment.completed,...
            ~isempty(experiment.err)};
    end
%% Callbacks
 
    function xlim_changed(~,eventdata)
        % Find fit line and redraw with more points
        ax_changed = eventdata.AffectedObject;
        uifitpeaks_lines = findobj(ax_changed,'tag','uifitpeaks');
        nlines = length(uifitpeaks_lines);
        xlim = ax_changed.XLim;
        for i = nlines:-1:1
            if isa(uifitpeaks_lines(i).UserData,'cfit')
                x = linspace(xlim(1),xlim(2),length(uifitpeaks_lines(i).XData));
                uifitpeaks_lines(i).XData = x;
                uifitpeaks_lines(i).YData = uifitpeaks_lines(i).UserData(x);
                return
            end
        end
    end
    function selector_click_callback(hObj,eventdata)
        if isempty(eventdata.Indices) || eventdata.Indices(2)~=10
            return
        end
        exp_ind = hObj.Data{eventdata.Indices(1),3};
        err = sites(site_index).experiments(exp_ind).err;
        if ~isempty(err)
            errmsg = getReport(err,'extended','hyperlinks','off');
            errmsg = strrep(errmsg,[newline newline],newline);
            f = figure('name',sprintf('Error (site: %i, exp: %i)',site_index,exp_ind),...
                'numbertitle','off','menubar','none','toolbar','none');
            uicontrol(f,'units','normalized','position',[0,0,1,1],'style','edit',...
                'string',errmsg,'max',Inf,...
                'HorizontalAlignment','left');
        end
    end
    function selector_edit_callback(hObj,eventdata)
        switch eventdata.Indices(2)
            case 1 % Display
                if busy
                    hObj.Data{eventdata.Indices(1),1} = ~hObj.Data{eventdata.Indices(1),1};
                    return
                end
                exp_ind = hObj.Data{eventdata.Indices(1),3};
                exp_type = hObj.UserData;
                mask = analysis.sites(site_index,exp_type).ignore==exp_ind;
                if any(mask) % Remove it
                    analysis.sites(site_index,exp_type).ignore(mask) = [];
                else % Add it
                    analysis.sites(site_index,exp_type).ignore(end+1) = exp_ind;
                end
                update_exp{exp_type}();
            case 6 % Redo Request
                state = hObj.Data{eventdata.Indices(1),6};
                % Do it for all others
                for i = 1:length([hObj.Data{:,5}])
                    if i == eventdata.Indices(1); continue; end
                    hObj.Data{i,6} = state;
                end
        end
    end
%% UIfitpeaks adaptor
    function save_state()
        if ~isstruct(ax(1).UserData) || ~isfield(ax(1).UserData,'uifitpeaks_enabled')
            analysis.sites(site_index,1).fit = [];
            analysis.sites(site_index,1).amplitudes = NaN;
            analysis.sites(site_index,1).locations = NaN;
            analysis.sites(site_index,1).etas = NaN;
            analysis.sites(site_index,1).widths = NaN;
            analysis.sites(site_index,1).background = NaN;
            analysis.sites(site_index,1).index = NaN;
            % "uses" and "redo" can stay untouched
            continue
        end
        fit_result = ax(1).UserData.pFit.UserData;
        new_data = true;
        analysis.sites(site_index,1).index = inds(site_index);
        if ~isempty(fit_result)
            fitcoeffs = coeffvalues(fit_result);
            if strcmpi(FitType,'voigt')
                nn = (length(fitcoeffs)-1)/4; % 4 degrees of freedom per peak for voigt; subtract background
            else
                nn = (length(fitcoeffs)-1)/3; % 3 degrees of freedom per peak; subtract background
            end
            analysis.sites(site_index,1).fit = fit_result;
            analysis.sites(site_index,1).amplitudes = fitcoeffs(1:nn);
            analysis.sites(site_index,1).locations = fitcoeffs(nn+1:2*nn);

            % Add widths and etas to structure
            switch lower(FitType)
                case 'gauss'
                    analysis.sites(site_index,1).widths = fitcoeffs(2*nn+1:3*nn)*2*sqrt(2*log(2));
                    analysis.sites(site_index,1).etas = zeros(1,nn);
                case 'lorentz'
                    analysis.sites(site_index,1).widths = fitcoeffs(2*nn+1:3*nn);
                    analysis.sites(site_index,1).etas = ones(1,nn);
                case 'voigt'
                    analysis.sites(site_index,1).widths = fitcoeffs(2*nn+1:3*nn);
                    analysis.sites(site_index,1).etas = fitcoeffs(3*nn+2:4*nn+1);
                otherwise
                    error(sprintf('Unsupported FitType %s',FitType));
            end

            analysis.sites(site_index,1).background = fitcoeffs(3*nn+1);
        else
            analysis.sites(site_index,1).fit = [];
            analysis.sites(site_index,1).amplitudes = [];
            analysis.sites(site_index,1).locations = [];
            analysis.sites(site_index,1).etas = [];
            analysis.sites(site_index,1).widths = [];
            analysis.sites(site_index,1).background = [];
        end
    end
    function attach_uifitpeaks(ax,init,varargin)
        % Wrapper to attach uifitpeaks
        % Let uifitpeaks update keyboard fcn, but then wrap that fcn again
        if any(isnan(init.locations))
            uifitpeaks(ax,'fittype',FitType,varargin{:});
        else
            uifitpeaks(ax,'fittype',FitType,'init',init,varargin{:});
        end
        if fig.UserData.uifitpeaks_count == 1 % Only set on first creation
            fig.UserData.uifitpeaks_keypress_callback = get(fig,'keypressfcn');
            set([fig, selector],'keypressfcn',@keypress_wrapper);
        end
    end
    function keypress_wrapper(hObj,eventdata)
        % uifitpeaks doesn't use alt, so we will distinguish with that
        if length(eventdata.Modifier)==1 && strcmp(eventdata.Modifier{1},'alt')
            cycleSite(hObj,eventdata);
        else
            hObj.UserData.uifitpeaks_keypress_callback(hObj,eventdata);
        end
    end

end
