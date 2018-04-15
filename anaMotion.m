function anaMotion

% anaMotion: analyse data collected with runMotion
% Test types:
%   testType = 10: raised Gabor until response, for calibration
%   testType = 12: raised Gabor ending in trigger pulse, for calibration
%   testType = 13: fading laterally-shifted Gabor
%   testType = 14: pulsed laterally-shifted Gabor
%   testType = 15: paired Gabors with variable orientation
%   testType = 16: flash, then fading laterally-shifted Gabor
%		testType = 17: moving edge
%   testType = 18: moving grating
%   testType = 19: moving bar defined by edge asynchrony
%		testType = 20: moving bar defined by width
%		testType = 21: pulsed laterally-shifted bar
%   testType = 22: fading Gabor with variable orientation
%   testType = 23: fading laterally-shifted Gabor with moving bar
%   testType = 24: flash then Gabor with variable orientation
%   testType = 25: fading tilted Gabor with moving bar
% Copyright:
%   Alan W Freeman

	% Specify the data source
	%	m.file = 'MotionFast5.lab';
	
	% Specify data to be selected
	m.select.Date = {'12-May-2017', '19-May-2017'}; % KN
	m.select.Radius = .2; % bar width (deg)
	%	m.select.Speed = 3; % speed
	%	m.select.Subject = {'HH', 'KT', 'LL', 'MA', 'SL', 'SW'};
	m.select.Subject = 'KN';
	m.select.TestType = 20; % test stimulus type
	
	% Specify default metadata
	m.plot.format = 'basic'; % plot format
	%	m.plot.arg = {'clipping', 'off'}; % no clipping at the axes edges

	% Specify tasks to be performed
	switch 'psych'
		case 'check' % numerical check
			m.tasks = 'suc mean list';
			m.limit.Cont = [-.015, -.01];
			m.limit.Cont = [.01, .015];
			m.x = 'Cont'; m.y = 'Success';
		case 'chron' % chronometric function
			m.tasks = 'light mag bin mean plot stop set';
			%	m.tasks = 'light mag bin mean fitlm show anova';
			m.limit.TimeR = [-1, 1];
			m.x = 'Cont'; m.y = 'TimeR';
			m.axes = {'TestType', 'Subject', 'Speed'}; m.line = {'Light'};
			%	m.axes = {'TestType', 'Speed'}; m.line = {'Subject', 'Light'};
			%	m.line = {'Light'};
			g = [m.axes, m.line]; m.bin.group = m.axes; m.bin.bins = 5; % no. of bins
			m.mean.group = [g, 'Bin'];
			% m.plot.arg = {'o-', 'clipping', 'off'}; % symbols, no clipping
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			%	m.plot.subplot = [3, 3]; % subplot layout
			m.fitlm.arg = {'TimeR ~ Subject + Light + ContMag^2'};
			m.set.arg = {'xLim', [0, .06], 'xTick', [0, .03, .06], ...
				'yLim', [.2, .5], 'yTick', [.2, .35, .5]};
			%	m.set.arg = {'xLim', [0, .17], 'xTick', [0, .05, .1, .15], ...
			%	'yLim', [.32, .48], 'yTick', [.35, .4, .45]};
			%	m.set.arg = {'xLim', [0, .6], 'xTick', [0, .3, .6], ...
			%	'yLim', [.375, .5], 'yTick', [.4, .45, .5]};
		case 'chronMean' % chronometric function mean plus confidence intervals
			m.tasks = 'light mag bin mean save meanC bar plot set err';
			m.limit.TimeR = [-1, 1];
			m.x = 'Cont'; m.y = 'TimeR';
			m.bin.group = {'TestType', 'Speed'}; m.bin.bins = 5;
			%	m.bin.group = {'TestType', 'Speed', 'Subject', 'Light'};
			m.mean.group = {'TestType', 'Speed', 'Subject', 'Light', 'Bin'};
			m.save.name = 'Bar chronometric 30';
			m.meanC.group = {'TestType', 'Speed', 'Light'};
			m.bar.group = {'TestType', 'Speed'};
			m.axes = {'TestType'};
			m.plot.fun = @(d, m)bar(d.(m.y));
			m.set.arg = {'yLim', [.1, .5], 'yTick', [.1, .3, .5]};
		case 'collect' % collect latency difference data into one plot: use in anaOr
			m.tasks = 'join plot';
			m.join.source = {'Low', 'High'};
			m.x = 'Speed'; m.y = 'Lat';
			m.line = {'Source'};
		case 'cor' % correlation check
			m.tasks = 'cor';
			m.select.ContC = 1; % select trials with bar present
		case 'lat' % latency difference versus speed
			m.tasks = 'light mag bin mean lat change conf plot set add stop save';
			%	m.tasks = 'light mag bin mean ext interp change fitlm show anova';
			m.limit.TimeR = [-1, 1];
			m.x = 'Cont'; m.y = 'TimeR';
			g = {'TestType', 'Speed', 'Subject'};
			%	g = {'TestType', 'Speed', 'Light'};
			m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Light', 'Bin'];
			m.lat.group = g;
			m.change.arg = {'x', 'Speed', 'y', 'LatDif'};
			m.conf.group = {'TestType', 'Speed'};
			m.axes = {'TestType'}; m.line = {};
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			%	m.plot.subplot = [3, 3]; % subplot layout
			m.set.arg = ...
				{'xScale', 'log', 'xLim', [.8, 32], 'xTick', [1, 3, 10, 30], ...
				'yLim', [-.03, .13], 'yTick', [0, .05, .1]};
			m.add.func = ...
				@(d, m)errorbar(d.(m.x), d.(m.y), d.Conf(:, 1), d.Conf(:, 2));
			m.save.name = 'Low';
		case 'psych' % psychometric function
			m.tasks = 'suc light mag bin mean plot stop set';
			%	m.tasks = 'suc light mag bin mean fitlm show anova';
			m.x = 'Cont'; m.y = 'Success';
			m.axes = {'TestType', 'Subject', 'Speed'}; m.line = {'Light'};
			%	m.axes = {'TestType', 'Speed'}; m.line = {'Subject', 'Light'};
			%	m.line = {'Light'};
			g = [m.axes, m.line]; m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Bin'];
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			%	m.plot.subplot = [3, 3]; % subplot layout
			m.fitlm.arg = {'Success ~ Subject + Light + ContMag^4'};
			m.set.arg = {'yLim', [.45, 1], 'yTick', [.5, .75, 1], ...
			'xLim', [0, .055], 'xTick', [0, .025, .05]};
			%	'xLim', [0, .13], 'xTick', [0, .05, .1]};
			%	'xLim', [0, .6], 'xTick', [0, .3, .6]};
		case 'psychMean' % psychometric function mean plus confidence intervals
			m.tasks = 'suc light mag bin mean meanP bar plot set err';
			m.x = 'Cont'; m.y = 'Success';
			m.bin.group = {'TestType', 'Speed'}; m.bin.bins = 5;
			%	m.bin.group = {'TestType', 'Speed', 'Subject', 'Light'};
			m.mean.group = {'TestType', 'Speed', 'Subject', 'Light', 'Bin'};
			m.meanP.group = {'TestType', 'Speed', 'Light'};
			m.bar.group = {'TestType', 'Speed'};
			m.axes = {'TestType'};
			m.plot.fun = @(d, m)bar(d.(m.y));
			m.set.arg = {'yLim', [.5, 1], 'yTick', [.5, .75, 1]};
		case 'psychWhole' % psychometric function
			m.tasks = 'suc bin mean plot stop set';
			% Set grouping variables
			m.axes = {'TestType', 'Subject'}; m.line = {'ContC'};
			g = [m.axes, m.line]; m.bin.group = g; m.mean.group = [g, 'Bin'];
			% Set other metadata
			% m.bin.bins = 5;
			m.plot.subplot = [3, 3]; % subplot layout
			m.suc.cor = 0; % 1 for correct response, 0 for choice of lighter side
			% m.set.arg = {'xLim', [0, .05], 'xTick', [0, .025, .05], ...
			m.set.arg = {'yLim', [.4, 1], 'yTick', [.5, .75, 1]};
			m.x = 'Cont'; m.y = 'Success';
		case 'sim' % simulate model
			m.tasks = 'par sim plot set';
			% Set grouping variables
			m.axes = {'Cont', 'Dir'}; m.line = {'Neuron'};
			m.sim.group = m.axes;
			% Set other metadata
			m.plot.colourOrder = [0, 1, 0; 0, 0, 0; 0, 0, 1; 1, 0, 0];
				% green, black, blue, red
			m.set.arg = {'xLim', [-.05, .15], 'xTick', [0, .1]};
			% m.sim.type = 'conv'; % type of impulse response
			m.x = 'Time'; m.y = 'Resp';
		case 'scat' % scatter plot of reaction times
			m.tasks = 'light plot set';
			% Set grouping variables
			m.axes = {'TestType', 'Subject', 'Speed'}; m.line = {'Light'};
			% Set other metadata
			m.limit.TimeR = [-1, 1];
			m.plot.arg = {'o', 'clipping', 'off'};
			%	m.plot.subplot = [3, 3]; % subplot layout
			% m.set.arg = {'xLim', [-.06, .06], 'xTick', [-.06, 0, .06], ...
			m.set.arg = {'yLim', [0, 1], 'yTick', [0, .5, 1]};
			m.x = 'Cont'; m.y = 'TimeR';
		case 'speedC' % contrast sensitivity versus speed for chronometric function
			m.tasks = 'light mag bin mean ext interp change plot set';
			m.tasks = 'light mag bin mean ext interp change fitlm show anova';
			m.limit.TimeR = [-1, 1];
			m.x = 'Cont'; m.y = 'TimeR';
			g = {'TestType', 'Speed', 'Subject', 'Light'};
			%	g = {'TestType', 'Speed', 'Light'};
			m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Bin'];
			m.ext.group = {'TestType', 'Speed', 'Subject'}; % m.ext.group = 'Speed';
			m.interp.group = g;
			m.change.arg = {'x', 'Speed', 'y', 'Sens'};
			m.axes = {'TestType', 'Subject'}; m.axes = {'TestType'};
			m.line = {'Subject', 'Light'}; m.line = {'Light'};
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			%	m.plot.subplot = [3, 3]; % subplot layout
			m.set.arg = {'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
				'yLim', [0, 70], 'yTick', [0, 35, 70]};
			m.fitlm.arg = {'Sens ~ Subject + Light*Speed + Speed^2'};
			%	m.fitlm.arg = {'Sens ~ Subject + Light + Speed^2'};
		case 'speedCFast' % cont. sens. vs. speed for psych. fn. incl. 30 deg/s
			m.tasks = 'light mag bin mean ext interp norm change plot set';
			m.limit.TimeR = [-1, 1];
			m.x = 'Cont'; m.y = 'TimeR';
			g = {'TestType', 'Radius', 'Speed', 'Light'};
			m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Bin'];
			m.ext.group = {'Radius', 'Speed'};
			m.interp.group = g;
			m.norm.group = 'Radius';
			m.change.arg = {'x', 'Speed', 'y', 'Sens'};
			m.axes = {'TestType'}; m.line = {'Light'};
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			m.set.arg = { ...
				'xScale', 'log', 'xLim', [1, 30], 'xTick', [1, 3, 10, 30], ...
				'yScale', 'log', 'yLim', [.05, 2], 'yTick', [.1, .3, 1]};
		case 'speedP' % contrast sensitivity versus speed for psychometric function
			% m.tasks = 'suc light mag bin mean fitnlm stop plot set pred add';
			m.tasks = 'suc light mag bin mean interp change plot set';
			m.tasks = 'suc light mag bin mean interp change fitlm show anova';
			m.x = 'Cont'; m.y = 'Success';
			g = {'TestType', 'Speed', 'Subject', 'Light'};
			%	g = {'TestType', 'Speed', 'Light'};
			m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Bin'];
			m.interp.group = g;
			% m.fitnlm.group = g;
			% m.fitnlm.x = 'ContMag'; m.fitnlm.y = 'Success';
			% m.fitnlm.funFun = @model; m.fitnlm.b0 = .05;
			m.change.arg = {'x', 'Speed', 'y', 'Sens'};
			m.axes = {'TestType', 'Subject'}; m.axes = {'TestType'};
			m.line = {'Subject', 'Light'}; m.line = {'Light'};
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			%	m.plot.subplot = [3, 3]; % subplot layout
			m.set.arg = {'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
				'yLim', [0, 150], 'yTick', [0, 75, 150]};
			m.fitlm.arg = {'Sens ~ Subject + Light*Speed + Speed^2'};
			%	m.fitlm.arg = {'Sens ~ Subject + Light*Speed'};
			m.pred.group = g;
		case 'speedPFast' % cont. sens. vs. speed for psych. fn. incl. 30 deg/s
			m.tasks = 'suc light mag bin mean interp norm change plot set';
			m.x = 'Cont'; m.y = 'Success';
			g = {'TestType', 'Radius', 'Speed', 'Light'};
			m.bin.group = g; m.bin.bins = 5;
			m.mean.group = [g, 'Bin'];
			m.interp.group = g;
			m.norm.group = 'Radius';
			m.change.arg = {'x', 'Speed', 'y', 'Sens'};
			m.axes = {'TestType'}; m.line = {'Light'};
			m.plot.colourOrder = [0, 0, 1; 1, 0, 0];
			m.set.arg = { ...
				'xScale', 'log', 'xLim', [1, 30], 'xTick', [1, 3, 10, 30], ...
				'yScale', 'log', 'yLim', [.05, 2], 'yTick', [.1, .3, 1]};
		case 'summary' % summarise data table
			m.tasks = 'summary';
			m.summary.exc = {'Date', 'Dir', 'AltC', 'Resp2'};
			m.summary.lim = [2, 20];
	end
	
	% Define handles for task functions in this file
	m.anova.fun = @doAnova;
	m.bar.fun = @doBar;
	m.change.fun = @doChange;
	m.conf.fun = @doConf;
	m.cor.fun = @doCor;
	m.err.fun = @doErr;
	m.ext.fun = @doExt;
	m.fitlm.fun = @doFitlm;
	m.fitnlm.fun = @doFitnlm;
	m.interp.fun = @doInterp;
	m.join.fun = @doJoin;
	m.light.fun = @doLight;
	m.mag.fun = @doMag;
	m.meanC.fun = @doMeanC;
	m.meanP.fun = @doMeanP;
	m.lat.fun = @doLat;
	m.norm.fun = @doNorm;
	m.par.fun = @doPar;
	m.save.fun = @doSave;
	m.show.fun = @doShow;
	m.sig.fun = @doSig;
	m.sim.fun = @doSim;
	m.suc.fun = @doSuc;

	% Compile the data table
	[d, m] = readDataTable(m); % generate the data table
	d = fixData(d); % correct errors in the data table
	
	% Analyse
	d = analyse(d, m);
	switch 0 % debug?
		case 1 % debug
			disp(d(1, :)); disp(d(end, :));
	end

function hErrorbar = addErr(h, y, e)

% Add error bars to an existing bar plot
%
% Inputs:
%   h = axes handle for existing plot
%   y = bar heights
%		e = error bar lengths
% Output:
%   hErrorbar = handle for error bars

	% Initialise
	b = h.Children; % handles of bar objects
	[rows, cols] = size(y); % size of bar height matrix
	hold on; % retain existing plotted data

	% Separate lower and upper errors
	if ndims(e) == ndims(y) + 1 % errors are potentially asymmetric
		eL = e(:, :, 1); eU = e(:, :, 2);
	elseif isvector(y) ~= isvector(e) % differing numbers of heights and errors
		eL = e(:, 1); eU = e(:, 2);
	else
		eL = e; eU = e;
	end

	if rows > 1
		hErrorbar = zeros(1,cols);
		for col = 1:cols

			% Extract the x location data needed for the errorbar plots
			x = b(col).XData - [b(col).XOffset]; % why minus?

			% Use the mean x values to call the standard errorbar function
			hErrorbar(col) = errorbar(mean(x, 1), y(:, col), ...
				eL(:, col), eU(:, col), '.k');
			set(hErrorbar(col), 'marker', 'none')

		end
	else
		x = b.XData - [b.XOffset];
		hErrorbar = errorbar(mean(x, 1), y, eL, eU, '.k');
		set(hErrorbar, 'marker', 'none')
	end
	
	% Finalise
	hold off;

function [d, m] = doAnova(d, m)

% Perform and display analyses of variance on the models

	for i = 1: length(m.model)
		disp(anova(m.model{i}));
	end

function [d, m] = doBar(d, m)

% Prepare bar plot: assume two rows with dark followed by light

	y = d.(m.y); c = d.Conf; % response and confidence interval
	d = d(1, :); % keep only first line
	d.(m.y) = y';
	d.Conf = reshape(c, [1, size(c)]);
	
function [d, m] = doChange(d, m)

% Change variables

	% Change variables in m.change.arg
	if isfield(m.change, 'arg') % m.change.arg is defined
		s = m.change.arg; % cell array of parameter/value pairs
		for i = 2: 2: length(s) % loop over pairs
			m.(s{i - 1}) = s{i};
		end
	end

function [d, m] = doConf(d, m)

% Calculate mean and confidence interval assuming Gaussian probability density

	muX = mean(d.(m.x)); % mean of x variable
	[muY, ~, muConf] = normfit(d.(m.y)); % mean and c.i. of response variable
	d = d(1, :); % keep only first row
	d.(m.x) = muX; d.(m.y) = muY; % means
	d.Conf = muConf' - muY; % confidence interval lengths
	
function [d, m] = doCor(d, m)

% Calculate correlation between conditioning and test alternatives

	[r, p] = corrcoef(d.Dir, d.AltC);
	fprintf('Correlation coefficient, p: %.3g, %.3g\n', r(1, 2), p(1, 2));

function [d, m] = doExt(d, m)

% Calculate minimum and maximum

	ext = [min(d.(m.y)), max(d.(m.y))]; % extremes
	d.Ext = repmat(ext, [size(d, 1), 1]); % store
	
function [d, m] = doFitlm(d, m)

% Fit linear model

	% Initiliase
	if isfield(m.fitlm, 'arg') % m.fitlm.arg specified
		arg = m.fitlm.arg; % user-specificed arguments
	else
		arg = {}; % default
	end
	
	% Remove unused categories
	for name = d.Properties.VarNames % loop over variables in d
		nameC = char(name); % name as string
		val = d.(nameC); % values of variable
		if iscategorical(val) % if this variable is categorical
			d.(nameC) = removecats(val); % remove unused categories
		end
	end
	
	% Fit
	m.model{m.group} = fitlm(d, arg{:});

function [d, m] = doErr(d, m)

% Add error bars to a bar plot

	h = m.handle.axes(m.group); % handle of current axes
	axes(h); % make axes current
	addErr(gca, d.(m.y), d.Conf);
	
function [d, m] = doFitnlm(d, m)

% Fit nonlinear model

	x = double(d(:, m.fitnlm.x)); % predictor variables
	y = d.(m.fitnlm.y); % response variable
	fun = m.fitnlm.funFun; % handle of model function
	b0 = m.fitnlm.b0; % initial estimate of model coefficients
	if isfield(m.fitnlm, 'arg') % m.fitnlm.arg specified
		arg = m.fitnlm.arg; % name/value pairs
	else
		arg = {}; % default is no pairs
	end
	m.model{m.group} = fitnlm(x, y, fun, b0, arg{:});

function [d, m] = doInterp(d, m)

% Interpolate psychometric or chronometric function

	% Initialise
	if ismember('Ext', d.Properties.VarNames) % chronometric function
		yq = mean(d.Ext(1, :)); % criterion y value
	else % psychometric function
		yq = .75;
	end
	
	% Interpolate
	d = sortrows(d, m.y); % sort y into ascending order
	y = d.(m.y); % obtain y
	i = find(y <= yq, 1, 'last'); % index of y immediately below yq
	j = find(y > yq, 1); % index of y immediately above yq
	d = d([i, j], :); % keep only those two rows
	x = d.(m.x); y = d.(m.y); % obtain x, y
	if size(x, 1) > 1 % can interpolate
		xq = interp1(y, x, yq); % interpolated contrast
	else
		xq = x; % use only value
	end
	s = 1 / xq; % contrast sensitivity
	d = d(1, :); % keep only first line
	d.Sens = s; % store
	d = describe(d, 'Sens', 'Contrast sensitivity'); % describe new variable

function [d, m] = doJoin(d, m)

% Collect latency difference data into one table

	s = m.join.source; % data sources
	for i = 1: length(s) % loop over sources
		name = s{i}; % current source
		load(name, 'd'); % current file
		name = repmat({name}, [size(d, 1), 1]); % one row for each row of d
		d.Source = categorical(name); % label data
		if isa(d, 'dataset') % d is a dataset
			d = dataset2table(d); % convert it to a table
		end
		if i == 1 % first pass
			dC = d; % initialise
		else
			dC = outerjoin(dC, d, 'mergeKeys', true); % join with previous sources
		end
	end
	d = dC; % store
	
function [d, m] = doLat(d, m)

% Calculate difference between dark and light reaction times

	i = max(d.Bin); % maximum bin number
	i = ceil(.5 * i);% middle bin
	d = d(d.Bin == i, :); % select rows containing middle bin
	l = d.TimeR(d.Light == 0) - d.TimeR(d.Light == 1); % latency (s)
	d = d(1, :); % keep only first row
	d.LatDif = l; % store
	d = describe(d, 'LatDif', 'Latency difference (s)'); % describe new variable
	
function [d, m] = doLight(d, m)

% Indicate whether stimulus is light or dark

	d.Light = nan(length(d), 1); % define variable
	i = ismember(d.TestType, [13, 17, 20, 21, 22, 23, 25]);
		% determined by contrast
	d.Light(i) = d.Cont(i) >= 0;
	i = d.TestType == 19; % light determined by asynchrony
	d.Light(i) = d.DurAsync(i) >= 0;
	
function [d, m] = doMag(d, m)

% Calculate magnitude of x-value

	switch m.x
		case 'Cont'
			d.ContMag = abs(d.Cont); % contrast magnitude
			d = describe(d, 'ContMag', 'Contrast magnitude'); % describe new variable
			m.x = 'ContMag'; % reset the independent variable
		case 'DurAsync'
			d.DurAsyncMag = abs(d.DurAsync);
			d = describe(d, 'DurAsyncMag', 'Asychrony magnitude (s)');
			m.x = 'DurAsyncMag';
	end

function [d, m] = doMeanC(d, m)

% Calculate mean response at middle bin of pchronometric function

	bins = max(d.Bin); bin = ceil(.5 * bins); % middle bin
	d = d(d.Bin == bin, :); % select that bin
	y = d.(m.y); % response
	[muY, ~, muConf] = normfit(y); % mean and c.i. of response variable
	d = d(1, :); % keep only first row
	d.(m.y) = muY; % means
	d.Conf = muConf' - muY; % confidence interval lengths

function [d, m] = doMeanP(d, m)

% Calculate mean response at middle bin of psychometric function

	%	y = d.(m.y); % response
	bins = max(d.Bin); bin = ceil(.5 * bins); % middle bin
	d = d(d.Bin == bin, :); % select that bin
	%	y = y(d.Bin == bin);
	%	[muY, ~, muConf] = normfit(y); % mean and c.i. of response variable
	[muY, muConf] = binofit(sum(d.Sum), sum(d.Rows));
	d = d(1, :); % keep only first row
	d.(m.y) = muY; % means
	d.Conf = muConf - muY; % confidence interval lengths
	
function [d, m] = doNorm(d, m)

% Normalise sensitivity by value for light bars moving at 3 deg/s

	i = d.Speed == 3 & d.Light == 1; % light bars moving at 3 deg/s
	s = d.Sens(i); % sensitivity
	d.Sens = d.Sens / s; % normalise
	i = d.Speed == 3 & d.Radius == .2; % fast monitor data at 3 deg/s
	d = d(~ i, :); % remove these data

function [d, m] = doPar(~, m)

% Define simulation parameters

	% Set values
	Speed = 10; % stimulus speed (deg/s)
	Stages = 4; % number of subcortical stages
	Separation = .1; % separation between on- and off-channels (deg)
	TimeConstOff = .009; % off-neuron time constant (s)
	TimeConstOn = .011; % on-neuron time constant (s)

	% Store the values in a data table
	Dir = [-1; 1; -1; 1]; % stimulus direction: -1 for leftward, 1 for rightward
	Cont = [1; 1; -1; -1]; % contrast
	rows = length(Dir); % number of rows in data table
	Speed = repmat(Speed, [rows, 1]);
	Separation = repmat(Separation, [rows, 1]);
	Stages = repmat(Stages, [rows, 1]);
	TimeConstOff = repmat(TimeConstOff, [rows, 1]);
	TimeConstOn = repmat(TimeConstOn, [rows, 1]);
	d = dataset(Cont, Dir, Speed, Separation, Stages, TimeConstOff, TimeConstOn);

function [d, m] = doSave(d, m)

% Save the data table

	save(m.save.name, 'd');

function [d, m] = doShow(d, m)

% Show models

	for i = 1: length(m.model)
		disp(m.model{i});
	end

function [d, m] = doSig(d, m)

% Test statistical significance

	% Initialise
	if isfield(m.sig, 'arg') % arguments specified?
		arg = m.sig.arg{:}; % yes
	else
		arg = 'linear'; % no, use the default
	end
	
	% Stop fitlm from using empty categories
	for name = d.Properties.VarNames % loop over variables in d
		nameC = char(name); % name as string
		val = d.(nameC); % values of variable
		if iscategorical(val) % if this variable is categorical
			d.(nameC) = removecats(val); % remove empty categories
		end
	end
	
	% Fit a linear model and do the stats
	m.model{m.group} = fitlm(d, arg); % fit a linear model
	disp(m.model{m.group}); % display the model
	anova(m.model{m.group}) % display the analysis of variance

function [d, m] = doSim(d, m)

% Simulate model

	% Initialise
	sig = (.4 / sqrt(2)) / d.Speed;
		% s.d. of correlation between stimulus and receptive field (s)
	stages = d.Stages; % number of stages
	vel = d.Dir * d.Speed; % velocity (deg/s)
	ts = 101; % number of time samples
	t = linspace(-.15, .25, ts)'; % time (s)
	t0 = .5 * d.Separation / vel; % time from centre to subfield (s)
	
	% Calculate responses
	tau = d.TimeConstOff; % time constant (s)
	pOff = - d.Cont * sim(t - t0, sig, stages, tau, m);
	tau = d.TimeConstOn; % time constant (s)
	pOn = d.Cont * sim(t + t0, sig, stages, tau, m); % on
	pCort = pOn + pOff; % cortical
	pInt = max(0, pCort); % impulse rate
	pInt = .1 * cumsum(pInt); % integral of cortical potential
	
	% Store the values in the data table
	n = 4; % number of response types
	d = repmat(d, [n * ts, 1]); % replicate parameters
	off = repmat({'Off'}, [ts, 1]); on = repmat({'On'}, [ts, 1]); % neuron types
	cort = repmat({'Cort'}, [ts, 1]); int = repmat({'Int'}, [ts, 1]);
	d.Neuron = categorical([off; on; cort; int]); % store the values in d
	d.Time = repmat(t, [n, 1]);
	d.Resp = [pOff; pOn; pCort; pInt];

function [d, m] = doSuc(d, m)

% Calculate success variable

	% Default case
	d.Success = d.Resp2 == d.Dir; % success is correct response
	d = describe(d, 'Success', 'Proportion correct');
	if isfield(m.suc, 'cor') && m.suc.cor == 0 % success is choosing lighter side
		i = d.Cont < 0; % rows on which contrast is negative
		d.Success(i) = ~ d.Success(i); % negate success for those rows
		d = describe(d, 'Success', 'Proportion lighter');
	end

function d = fixData(d)
		
% Correct errors and omissions in the data table
	
	% Initialise
	d.Date = categorical(d.Date); % categorical
	d.ExpFunction = categorical(d.ExpFunction);
	d.Subject = categorical(d.Subject);
	d.Time = categorical(d.Time);
	
	% Remove bad rows
	d = d(d.ContTest ~= .7, :); % contTest incorrectly set
	d = d(d.Date ~= '05-May-2017', :); % RTBox incorrectly initalised
	i = d.Date == '15-May-2017' & d.Subject == 'LC'; d = d(~ i, :);
		 % didn't reach threshold
	
	% Remove trials with no response
	t = d.Resp2 <= 0; % trials with no response
	d = d(~ t, :); % remove trials
	
	% Add variable descriptions
	d = describe(d, 'Cont', 'Contrast');
	d = describe(d, 'Speed', 'Stimulus speed (deg/s)');
	d = describe(d, 'TimeR', 'Reaction time (s)');

function [d, m] = readDataTable(m)

% Define input and output files, and read data

	% Specify files
	if isfield(m, 'file') % single input file
		in = m.file;
		out = ''; % a mat file would be a nuisance, so don't create it
	else % all lab files
		in = {
			'Motion1.lab', 'Motion2.lab', 'Motion3.lab', 'Motion4.lab', ...
			'Motion5.lab', 'Motion6.lab', 'Motion7.lab', 'Motion8.lab', ...
			'Motion9.lab', 'Motion10.lab', 'Motion11.lab', 'Motion12.lab', ...
			'Motion13.lab', 'Motion14.lab', 'Motion15.lab', 'Motion16.lab', ...
			'Motion17.lab', 'MotionFast1.lab', 'MotionFast2.lab', ...
			'MotionFast3.lab', 'MotionFast4.lab', 'MotionFast5.lab'};
		out = 'Motion.mat';
	end
	
	% Read data
	d = readData(in, out);

function p = sim(t, sig, n, tau, m)

% Simulate time course of postsynaptic potential in cortical neuron

	% Initialise
	if isfield(m.sim, 'type')
		type = m.sim.type; % user-supplied value
	else
		type = 'gamma'; % default
	end
	
	% Calculate
	switch type
		
		case 'biphasic' % derivative of gamma density
			p = .5 * (1 + sign(t)) .* ...
				((n - 1 - t / tau) / (factorial(n - 1) * tau ^ n)) .* ...
				t .^ (n - 2) .* exp(-t / tau);
			
		case 'conv' % convolution of Gaussian and gamma densities
			
			% Calculate Gaussian density
			dt = t(2) - t(1); % time increment (s)
			tLim = 3 * [- sig, sig]; % time limits (s)
			ts = (tLim(2) - tLim(1)) / dt; % number of times
			ts = 2 * ceil(.5 * ts) + 1; % make it odd
			tC = linspace(tLim(1), tLim(2), ts); % times at which to compute (s)
			r = normpdf(tC, 0, sig); % density
			r = r / normpdf(0, 0, sig); % normalised density
			
			% Calculate convolution
			g = gampdf(t, n, tau); % gamma density
			p = conv(g, r, 'same'); % convolution
			
		case 'gamma' % gamma density
			p = gampdf(t, n, tau);
			
	end
	