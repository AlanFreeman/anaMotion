function [d, m] = analyse(d, m)

% analyse: sequentially analyse a data table
%
% Call: [d, m] = analyse(d, m)
% Inputs and outputs:
%   d = data table: dataset
%   m = metadata: structure
% Description:
%   Rows in the data table represent observations and columns represent
%   variables. Data are analysed in two sequences. First, the tasks listed
%   in m.tasks are performed on the table in sequence. Second, for each
%   task, data in the table are separated into groups and each group is
%   analysed in turn. Grouping variables are contained in m.task.group.
%   Plotting tasks differ in that there are two grouping variables: m.axes
%   indicates variables that change between axes and m.line indicates
%   variables that change between lines on an axes.
% Copyright:
%   Alan W Freeman

	% Select data
	d = select(d, m); % select rows in which variables with specified values
	d = limit(d, m); % select rows in which variables meet specified limits

	% Provide default grouping variables
	if ~ isfield(m, 'axes'), m.axes = {}; end % single axes
	if ~ isfield(m, 'line'), m.line = {}; end % single line per axes
	
	% Define handles for task functions in this file
	m.bin.fun = @doBin;
	m.add.fun = @doAdd; m.add.group = m.axes;
	m.list.fun = @doList;
	m.mean.fun = @doMean;
	m.pred.fun = @evalModel;
	m.set.fun = @doSet;
	
	% Perform tasks
	if ~ isempty(m.tasks)
		task = strtrim(m.tasks); % remove insignificant white space
	else
		return % nothing to do
	end
	task = strsplit(task); % split the string into a cell array of strings
	for i = 1: length(task) % loop over tasks
		taskC = task{i}; % current task
		if strcmp(taskC, 'stop') % stop processing
			break
		end
		if isempty(d) % data table is empty
			error('Data table input to task ''%s'' is empty', taskC);
		end
		if strcmp(taskC, 'plot') % plotting function
			[d, m] = plotGroups(d, m);
		elseif strcmp(taskC, 'summary') % summarise data table
			[d, m] = doSummary(d, m);
		elseif isfield(m, taskC) && isfield(m.(taskC), 'fun')
			% function handle exists
			[d, m] = anaGroups(d, m, taskC);
		else % error
			error('Task ''%s'' is undefined', taskC);
		end
	end

function [d, m] = anaGroups(d, m, task)

% Analyse data group by group
%
% Input:
%   task = task currently being executed: string

	% Find grouping variable(s) and function to execute
	if isfield(m.(task), 'group') % grouping variable(s) for this task
		var = m.(task).group;
	else % no grouping
		var = {};
	end
	fun = m.(task).fun; % handle of function to execute

	% Split d into groups
	g = group(d, var); % group data
	gs = length(g); % number of groups
	
	% Loop over groups
	d = cell(gs, 1); % allocate storage
	for j = 1: gs
		m.group = j; % provide group number
		[d{j}, m] = fun(g{j}, m); % analyse group
	end
	
	% Finalise
	d = vertcat(d{:}); % combine data tables over groups
	
function [d, m] = countTrials(d, m)

% Count rows in the data table

	rows = size(d, 1); % number of rows
	d = d(1, :); % keep only first row
	d.Rows = rows; % add variable

function [d, m] = doAdd(d, m)

% Add data to existing axes

	% Initialise
	persistent fun gLine;
	if m.group == 1 % first time through
		if isfield(m.add, 'func') % set plotting function
			fun = m.add.func; % user-provided function
		else
			if isfield(m.add, 'arg') % plot arguments
				arg = m.add.arg; % user-supplied arguments
			else
				arg = {}; % default
			end
			fun =  @(d, m)plot(d.(m.x), d.(m.y), arg{:}); % default
		end
		gLine = m.line; % variables that change between lines
	end

	% Loop over lines
	hC = m.handle.axes(m.group); % axes handle
	axes(hC); % make axes current
	aC = d; % current data table
	p = group(aC, gLine); % group into lines
	ps = length(p); % number of lines
	set(hC, 'nextPlot', 'add'); % add lines, don't replace
	if isfield(m, 'plot') && isfield(m.plot, 'colourOrder') % pre-R2014b
		colourOrder = m.plot.colourOrder; % line colours
	else
		colourOrder = get(hC, 'ColorOrder');
	end
	%{
	% set line colours in R2014b+
	if isfield(m, 'plot') && isfield(m.plot, 'colourOrder')
		hC.colorOrder = m.plot.colourOrder;
	end
	set(hC, 'colorOrderIndex', 1); % start at first colour
	%}
	%	hLine = findobj(hC, 'type', 'line'); % handles of existing lines
	%	set(hLine, 'lineStyle', 'none', 'marker', 'o');
		% reset existing line style to distinguish from model
	for i = 1: ps % loop over lines
		pC = p{i}; % current data to plot
		hP = fun(pC, m); % plot: pre-R2014b
		j = mod(i - 1, size(colourOrder, 1)) + 1; % set line color
		set(hP, 'color', colourOrder(j, :));
		% fun(pC, m); % plot: R2014b+
	end
	
function [d, m] = doBin(d, m)

% Calculate bin number, with (as far as possible) equal numbers of samples
% per bin

	if isfield(m.bin, 'bins') % set number of bins
		bins = m.bin.bins; % user-supplied
	else
		bins = 10; % default
	end
	ts = size(d, 1); % number of trials
	d = sortrows(d, m.x); % sort by stimulus value
	d.Bin = (1 + floor(bins * ((0: ts - 1) / ts)))'; % bin number of each trial

function [d, m] = doList(d, m)

% List data table

	if isfield(m.list, 'vars') % variables to list
		disp(d(:, m.list.vars));
	else % list whole data table
		disp(d);
	end

function [d, m] = doMean(d, m)

% Calculate mean stimulus and response

	stimMean = mean(d.(m.x)); % mean stimulus
	rows = size(d, 1); % number of rows
	respSum = sum(d.(m.y)); % summed response
	respMean = mean(d.(m.y)); % mean response
	d = d(1, :); % keep only first row
	d.(m.x) = stimMean;
	d.Rows = rows;
	d.Sum = respSum;
	d.(m.y) = respMean;

function [d, m] = doSet(d, m)

% Set properties of existing axes

	if isfield(m.set, 'arg')
		set(m.handle.axes, m.set.arg{:});
	end

function [d, m] = doSummary(d, m)

% List variables with specified numbers of unique values. Add a variable
% giving the numbers of rows with each set of unique values.

	% Exclude specified variables
	if isfield(m, 'summary') && isfield(m.summary, 'exc')
		var = d.Properties.VarNames; % names
		var = ~ ismember(var, m.summary.exc); % find variables that are staying
		d = d(:, var); % remove specified variables
	end
	
	% Find variables with specified numbers of unique values
	if isfield(m, 'summary') && isfield(m.summary, 'lim')
		% limits on number of unique values of a variable
		lim = m.summary.lim; % user-supplied value
	else
		lim = [2, 10]; % default
	end
	vs = size(d, 2); % number of variables
	var = false(vs, 1); % variables to keep
	for v = 1: vs % loop over variables
		dC = d(:, v); % current variable
		dC = unique(dC); % unique values for this variable
		rowsC = size(dC, 1); % number of unique rows
		if rowsC >= lim(1) && rowsC <= lim(2) % number lies within specified limits
			var(v) = true; % add this variable
		end
	end
	var = d.Properties.VarNames(var); % names
	d = d(:, var); % keep only those variables
	
	% Calculate number of rows
	m.count.fun = @countTrials;
	m.count.group = var;
	[d, m] = anaGroups(d, m, 'count');
	
	% List
	disp(d);

function [d, m] = evalModel(d, m)

% Evaluate model

	d.(m.y) = predict(m.model{m.group}, d.(m.x));

function g = group(d, col)
		
% Group data table into rows with unique parameters
%
% Input:
%   d = data table
%   col = cell array of names of parameters that are constant within a group
% Output:
%   g = d split into a cell array, with one member for each group

	if isempty(col) % no grouping
		g{1} = d;
	else
		u = unique(d(:, col)); % unique rows
		gs = size(u, 1); % number of groups
		g = cell(1, gs); % allocate storage
		for i = 1: gs % loop over groups
			r = ismember(d(:, col), u(i, :)); % find all rows in this group
			g{i} = d(r, :); % assign them to this group
		end
	end
	
function [h, i] = iniAxes(i, m)

% Initialise new axes
%
% Input:
%   i = axes index within complete data set
%   m = metadata structure
% Output:
%   h = axes handle
%   i = axes index within figure
%	Notes:
%		plot formats are 'basic', 'paper', 'poster', 'talk', or 'thesis'

	% Default input arguments
	if isfield(m, 'plot') && isfield(m.plot, 'format') % plot format
		form = m.plot.format;
	else % default
		form = 'basic';
	end
	if isfield(m, 'plot') && isfield(m.plot, 'subplot') % subplot layout
		layout = m.plot.subplot;
	else % default
		if strcmp(form, 'basic')
			layout = [2, 2];
		else
			layout = [1, 1];
		end
	end
	axesPerFig = prod(layout); % number of axes per figure
	
	% Create axes
	i = mod(i - 1, axesPerFig) + 1; % axes index on this figure
	if i == 1
		figure('windowStyle', 'docked'); % new figure
	end
	h = subplot(layout(1), layout(2), i);

	% Set formatting parameters
	switch form
		case 'paper'
			set(h, 'fontSize', 16); % axis label size (pt)
		case 'poster'
			set(h, 'fontSize', 28);
		case 'talk'
			set(h, 'fontSize', 16);
		case 'thesis'
			set(h, 'fontSize', 12);
	end

function labelAxes(d, gAxes)

% Label axes with values of axes variables

	vs = length(gAxes); % number of variables
	if vs == 0 % no label required
		return
	end
	s = cell(1, vs); % one string for each variable
	for i = 1: vs
		val = d.(gAxes{i})(1, :); % value of current axes variable
		s{i} = makeString(val); % convert to char
	end
	s = sprintf('%s, ', s{:}); % concatenate
	s = s(1: end - 2); % remove trailing punctuation
	title(s);

function labelAxis(d, a)

% Label x and y axes
%	Input:
%		d = data table
%		a = cell array of x- and y-variable names

	n = cell(2, 1); % x- and y-names
	des = d.Properties.VarDescription; % variable descriptions
	desC = ''; % current description
	for j = 1: 2 % x- then y-axes
		i = strcmp(a{j}, d.Properties.VarNames);  % variable index
		if ~ isempty(des) % there is at least one description
			desC = des{i}; % variable description
		end
		if isempty(desC) % no description supplied
			n{j} = a{j}; % use variable name
		else
			n{j} = desC; % use description
		end
	end
	xlabel(n{1}); ylabel(n{2});

function labelFigure(i, gAxes, gLine)

% Label figure: axes variables followed by line variables

	if i > 1 % figure is already labelled
		return
	end
	sAxes = []; sLine = []; % initialise
	if ~ isempty(gAxes) % there is at least one axes variable
		sAxes = sprintf('%s, ', gAxes{:}); % list variables
		sAxes = sAxes(1: end - 2); % remove trailing punctuation
	end
	if ~ isempty(gLine) % there is at least one line variable
		sLine = sprintf('%s, ', gLine{:}); % list variables
		sLine = sLine(1: end - 2); % remove trailing punctuation
	end
	if ~ isempty(sAxes) && ~ isempty(sLine) % there are both axes and line var.
		s = [sAxes '; ' sLine]; % concatenate
	else
		s = [sAxes, sLine]; % no need for separator
	end
	if ~ isempty(s) % insert label at top left
		annotation('textbox', [0 .9 .1 .1], 'string', s, 'lineStyle', 'none');
	end

function labelLine(d, h, gLine)

% Label line: add its display name to the line's lineseries object

	vs = length(gLine); % number of variables
	if vs == 0 % no label required
		return
	end
	s = cell(1, vs); % one string for each variable
	for i = 1: vs
		val = d.(gLine{i})(1, :); % value of current line variable
		s{i} = makeString(val); % convert to char
	end
	s = sprintf('%s, ', s{:}); % concatenate
	s = s(1: end - 2); % remove trailing punctuation
	set(h, 'displayName', s); % display name

function d = limit(d, m)
		
% Select rows in which specified variables fall within specified limits
%
% Input:
%		Variables are specified by fields in the metadata struture, m.limit.
%		Allowable values are specified as a closed range.

	% Is selection required?
	if isfield(m, 'limit') % yes
		
		% Initialise
		field = fieldnames(m.limit); % names of fields in m.limit
		fields = length(field); % number of fields
		rows = length(d); % number of rows
		row = zeros(rows, fields); % rows to select

		% For each field, find rows to select
		for i = 1: fields % loop over fields
			name = field{i}; % name of current field
			var = d.(name); % values of variable
			val = m.limit.(name); % value of field
			classC = class(var); % class of variable
			switch classC % check that the variable is double
				case 'double'
					rowC = var >= val(1) & var <= val(2);
				otherwise
					warning('Row selection: class %s not supported', classC);
			end
			row(:, i) = rowC; % store result for this field
		end

		% Select rows
		row = all(row, 2); % find rows that satisfy all tests
		d = d(row, :); % select rows
		
	end

function s = makeString(val)

% Convert value, numeric or not, to char

	if islogical(val) % logical value
		val = double(val); % convert to double
	end
	if isnumeric(val) % numeric value, including row vector
		s = sprintf('%g ', val); % convert to char
		s = s(1: end - 1); % remove trailing space
	else % char, cell, or categorical
		s = char(val); % convert to char if necessary
	end

function [d, m] = plotGroups(d, m)

% Plot data defined by grouping variables

	% Initialise
	gAxes = m.axes; % variables that change between axes
	gLine = m.line; % variables that change between lines
	if isfield(m, 'plot') && isfield(m.plot, 'fun') % set plotting function
		fun = m.plot.fun; % user-provided function
	else
		if isfield(m, 'plot') && isfield(m.plot, 'arg') % plot arguments
			arg = m.plot.arg; % user-supplied arguments
		else
			arg = {}; % default
		end
		fun =  @(d, m)plot(d.(m.x), d.(m.y), arg{:}); % default
	end

	% Loop over axes
	a = group(d, gAxes); % group into axes
	as = length(a); % number of axes
	h = zeros(1, as); % axes handles
	for j = 1: as % loop over axes

		% Loop over lines
		aC = a{j}; % current axes
		p = group(aC, gLine); % group into lines
		ps = length(p); % number of lines
		[h(j), i] = iniAxes(j, m); % initialise axes
		labelFigure(i, gAxes, gLine); % label figure with variable names
		labelAxes(aC, gAxes); % label axes with variable values
		if isfield(m, 'plot') && isfield(m.plot, 'colourOrder')
			set(h(j), 'colorOrder', m.plot.colourOrder); % set line colours
		end
		hold all; % add plots, don't replace
		for i = 1: ps % loop over lines
			
			% Plot a line
			pC = p{i}; % current data to plot
			s = func2str(fun); % convert function handle to string
			if strfind(s, 'contour') % function is contour
				[~, hP] = fun(pC, m);
			else % assume function is one, such as plot, with a single output
				hP = fun(pC, m);
			end
			labelLine(pC, hP, gLine); % label line by adding display name
			
		end
		if ~ isempty(gLine) % add legend
			m.handle.legend(j) = legend('location', 'best'); legend('boxOff');
		end
		labelAxis(aC, {m.x, m.y}); % label x- and y-axes
		hold off;
		
	end
	prettyPlot(h, m);
	m.handle.axes = h; % store axes handle(s) in m

function prettyPlot(h, m)

% Beautify an existing plot
%
% Input:
%   h = axes handle(s)
%		m = metadata structure
%	Notes:
%   plot formats are 'basic', 'paper', 'poster', 'talk', or 'thesis'

	% Default input arguments
	if isfield(m, 'plot') && isfield(m.plot, 'format')
		form = m.plot.format;
	else % default
		form = 'basic';
	end
	
	% Set basic plot parameters
	set(h, 'plotBoxAspectRatio', [1 1 1], 'box', 'off');
	set(h, 'tickLength', [.02 .02]);
	set(h, 'fontName', 'Helvetica');
	if strcmp(form, 'basic') % basic format
		return % don't do anything else
	end
	
	% Define advanced plot parameters
	pos = [5, 5, 5, 5]; % last two digits set axis lengths (cm)
	switch form
		case 'paper'
			wid = 1.5; % line width (pt)
		case 'poster'
			wid = 2;
		case 'talk'
			wid = 2;
		case 'thesis'
			wid = 1;
	end
	
	% Set advanced plot parameters
	set(h, 'units', 'centimeters', 'position', pos); % location and size
	set(h, 'lineWidth', wid); % line width (pt)
	% set(h, 'FontSize', 20); % axis title size (pt)
	% set(h, 'FontSize', 16); % axis label size (pt)
	hLine = findobj(h, 'type', 'line'); % handles of all lines
	set(hLine, 'lineWidth', wid); % set line width

function d = select(d, m)
		
% Select rows in which specified variables have specified values.
%
% Input:
%		Values are specified by fields in the metadata struture, m. Fields and
%		variables are matched by having the same name. Allowable values are:
%			char-based: categorical, cell or char
%			double: single values or a closed range

	% Is selection required?
	if isfield(m, 'select') % yes
		
		% Initialise
		field = fieldnames(m.select); % names of fields in m.select
		fields = length(field); % number of fields
		rows = length(d); % number of rows
		row = zeros(rows, fields); % rows to select

		% For each field, find rows to select
		for i = 1: fields % loop over fields
			name = field{i}; % name of current field
			val = m.select.(name); % value of field
			var = d.(name); % values of variable
			classC = class(var); % class of variable
			switch classC % selection depends on class
				case {'categorical', 'cell', 'char'}
					rowC = ismember(var, val);
				case 'double'
					%{
					if size(var(1, :), 2) == 1 && size(val, 2) > 1
						% variable is scalar but field is not: assume range
						rowC = var >= val(1) & var <= val(2);
					else
						rowC = ismember(var, val, 'rows');
					end
					%}
					rowC = ismember(var, val, 'rows');
				otherwise
					warning('Row selection: class %s not supported', classC);
			end
			row(:, i) = rowC; % store result for this field
		end

		% Select rows
		row = all(row, 2); % find rows that satisfy all tests
		d = d(row, :); % select rows
		
	end
