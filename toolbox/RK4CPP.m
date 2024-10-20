function [t, solution] = RK4CPP(treeModel, tmax, h)
    % RK4CPP Runs the RK4 simulation using the mex solver
    %
    % Inputs:
    % - treeModel: DendriteTreeModel object
    % - tmax: Maximum simulation time
    % - h: Time step size
    %
    % Outputs:
    % - t: Time vector
    % - solution: Simulated state variables over time

    % Validate inputs
    if ~isa(treeModel, 'DendriteTreeModel')
        error('First input must be a DendriteTreeModel object.');
    end
    if ~isscalar(tmax) || tmax <= 0
        error('tmax must be a positive scalar.');
    end
    if ~isscalar(h) || h <= 0
        error('h must be a positive scalar.');
    end

    % Prepare time vector
    t = 0:h:tmax;
    numTimeSteps = length(t);

    % Get initial conditions
    X0 = getInitialConditions(treeModel);

    % Prepare dendrite parameters
    numDendrites = treeModel.numDendrites;
    paramsArray = zeros(numDendrites, 7); % [alpha, b, Tau, TauR, NaX, gc, ID]

    for idx = 1:numDendrites
        dendrite = treeModel.dendrites(idx);
        paramsArray(idx, :) = [
            dendrite.params.alpha, 
            dendrite.params.b, 
            dendrite.params.Tau, 
            dendrite.params.TauR, 
            dendrite.params.NaX, 
            dendrite.params.gc, 
            dendrite.ID
        ];
    end

    % Prepare coupling matrix data
    [rows, cols, values] = findCouplingMatrix(treeModel);
    couplingData = [rows, cols, values]; % [numCouplings x 3]

    % Prepare StimuliMatrix and stimDendriteIDs
    [StimuliMatrix, stimDendriteIDs] = prepareStimuliMatrix(treeModel, numTimeSteps);

    % Prepare stimDendriteIDArray
    stimDendriteIDArray = stimDendriteIDs(:); % Column vector

    % Prepare dendriteIDArray
    dendriteIDArray = [treeModel.dendrites.ID]';
    
    % SOMA ONLY
    if treeModel.numDendrites == 1
        couplingData = [0, 0, 0];
    end

    % Run the MEX function
    [t, solution] = mexRK4Solver(paramsArray, couplingData, StimuliMatrix, stimDendriteIDArray, X0, h, numTimeSteps, dendriteIDArray);
end

%% Helper function to get initial conditions
function X0 = getInitialConditions(treeModel)
    numDendrites = treeModel.numDendrites;
    X0 = zeros(numDendrites * 2, 1); % [u0; v0] for each dendrite
    for i = 1:numDendrites
        X0(2 * i - 1) = treeModel.dendrites(i).InitialConditions(1);
        X0(2 * i)     = treeModel.dendrites(i).InitialConditions(2);
    end
end

%% Helper function to find coupling matrix entries
function [rows, cols, values] = findCouplingMatrix(treeModel)
    numDendrites = treeModel.numDendrites;
    rows = [];
    cols = [];
    values = [];

    for i = 1:numDendrites
        conn = treeModel.Connectivity(i);
        gc_i = treeModel.dendrites(i).params.gc;
        % Proximal connections
        for pid = conn.List_Proximal
            sourceIdx = i - 1; % Zero-based indexing
            targetIdx = treeModel.getDendriteIndex(pid) - 1; % Zero-based
            rows(end+1, 1) = sourceIdx;
            cols(end+1, 1) = targetIdx;
            values(end+1, 1) = gc_i;
        end
        % Distal connections
        for did = conn.List_Distal
            sourceIdx = i - 1; % Zero-based indexing
            targetIdx = treeModel.getDendriteIndex(did) - 1; % Zero-based
            rows(end+1, 1) = sourceIdx;
            cols(end+1, 1) = targetIdx;
            values(end+1, 1) = gc_i;
        end
    end
end

%% Helper function to prepare StimuliMatrix
function [StimuliMatrix, stimDendriteIDs] = prepareStimuliMatrix(treeModel, numTimeSteps)
    stimDendriteIDs = [];
    stimDendrites = findStimulatedDendrites(treeModel);
    numStimDendrites = length(stimDendrites);

    for sdIdx = 1:numStimDendrites
        stimDendriteIDs(end+1) = stimDendrites{sdIdx}.ID;
    end
    % Initialize StimuliMatrix as a flat matrix [numStimDendrites x numTimeSteps]
    StimuliMatrix = zeros(numStimDendrites, numTimeSteps);

    for i = 1:numStimDendrites
        StimSignal = treeModel.Stimuli{treeModel.getDendriteIndex(stimDendriteIDs(i))};
        % Assuming StimSignal is a vector of length numTimeSteps
        if length(StimSignal) ~= numTimeSteps
            error('Stimulus signal length does not match number of time steps.');
        end
        StimuliMatrix(i, :) = StimSignal;
    end
end

%% Helper function to find stimulated dendrites
function stimDendrites = findStimulatedDendrites(treeModel)
    stimDendrites = [];
    for i = 1:treeModel.numDendrites
        if ~isempty(treeModel.Stimuli{i})
            stimDendrites{end+1} = treeModel.dendrites(i);
        end
    end
end
