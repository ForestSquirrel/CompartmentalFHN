classdef DendriteTreeModel < handle
    properties
        dendrites       % Array of Dendrite objects
        numDendrites    % Number of dendrites
        Stimuli         % Cell array of stimuli
        Connectivity    % Struct array representing connections
        id_to_idx       % Mapping from ID to index
    end

    methods
        %% Constructor
        function obj = DendriteTreeModel()
            obj.dendrites = Dendrite.empty;
            obj.numDendrites = 0;
            obj.Stimuli = {};
            obj.Connectivity = struct('ID', {}, 'List_Proximal', {}, 'List_Distal', {});
            obj.id_to_idx = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
        end

        %% Add a Dendrite
        function addDendrite(obj, dendrite)
            if ~isa(dendrite, 'Dendrite')
                error('Input must be a Dendrite object.');
            end
            obj.numDendrites = obj.numDendrites + 1;
            obj.dendrites(obj.numDendrites) = dendrite;
            obj.Stimuli{obj.numDendrites} = [];
            obj.id_to_idx(dendrite.ID) = obj.numDendrites;
            % Initialize connectivity for this dendrite
            obj.Connectivity(obj.numDendrites).ID = dendrite.ID;
            obj.Connectivity(obj.numDendrites).List_Proximal = [];
            obj.Connectivity(obj.numDendrites).List_Distal = [];
        end

        %% Get Dendrite Index by ID
        function idx = getDendriteIndex(obj, ID)
            if isKey(obj.id_to_idx, ID)
                idx = obj.id_to_idx(ID);
            else
                error('Dendrite ID %d not found.', ID);
            end
        end

        %% Validate Dendrite ID
        function valid = isValidID(obj, ID)
            valid = isKey(obj.id_to_idx, ID);
        end

        %% Add Connection Between Dendrites
        function addConnection(obj, ID, List_Proximal, List_Distal)
            idx = obj.getDendriteIndex(ID);
            if ~isempty(List_Proximal)
                obj.Connectivity(idx).List_Proximal = unique([obj.Connectivity(idx).List_Proximal, List_Proximal]);
                % Add reciprocal connections
                for pid = List_Proximal
                    if obj.isValidID(pid)
                        pidx = obj.getDendriteIndex(pid);
                        obj.Connectivity(pidx).List_Distal = unique([obj.Connectivity(pidx).List_Distal, ID]);
                    else
                        error('Proximal Dendrite ID %d is invalid.', pid);
                    end
                end
            end
            if ~isempty(List_Distal)
                obj.Connectivity(idx).List_Distal = unique([obj.Connectivity(idx).List_Distal, List_Distal]);
                % Add reciprocal connections
                for did = List_Distal
                    if obj.isValidID(did)
                        didx = obj.getDendriteIndex(did);
                        obj.Connectivity(didx).List_Proximal = unique([obj.Connectivity(didx).List_Proximal, ID]);
                    else
                        error('Distal Dendrite ID %d is invalid.', did);
                    end
                end
            end
        end

        %% Add Connections via Connection String
        function addConnectionStr(obj, connStr)
            % Parse the input connection string and add connections to the model
            % according to the specified rules, allowing groups anywhere in the string.

            % Remove any spaces from the input string
            connStr = strrep(connStr, ' ', '');

            % Check if the string starts with '-'
            if startsWith(connStr, '-')
                error('Invalid connection string: Soma cannot have "-" to its left.');
            end

            % Split the string into tokens, respecting groups
            tokens = splitConnectionString(connStr);

            % Parse tokens into lists of IDs
            tokenIDs = cell(1, length(tokens));
            for i = 1:length(tokens)
                token = tokens{i};
                tokenIDs{i} = parseToken(token);
            end

            % Validate IDs
            allIDs = unique([tokenIDs{:}]);
            for i = 1:length(allIDs)
                if ~obj.isValidID(allIDs(i))
                    error('Invalid ID %d in connection string.', allIDs(i));
                end
            end

            % Now, for each pair of adjacent tokens, create connections
            for i = 1:length(tokenIDs)-1
                idsA = tokenIDs{i};
                idsB = tokenIDs{i+1};
                % Create connections from every ID in idsA to every ID in idsB
                for a = idsA
                    for b = idsB
                        obj.addConnection(a, [], b);
                    end
                end
            end

            %% Nested function to split the connection string
            function tokens = splitConnectionString(str)
                str = char(str);
                tokens = {};
                currentToken = '';
                insideGroup = false;
                idx = 1;
                while idx <= length(str)
                    ch = str(idx);
                    if ch == '['
                        insideGroup = true;
                        currentToken = [currentToken, ch];
                    elseif ch == ']'
                        insideGroup = false;
                        currentToken = [currentToken, ch];
                    elseif ch == '-' && ~insideGroup
                        % End of token
                        if ~isempty(currentToken)
                            tokens{end+1} = currentToken;
                            currentToken = '';
                        end
                    else
                        currentToken = [currentToken, ch];
                    end
                    idx = idx + 1;
                end
                if ~isempty(currentToken)
                    tokens{end+1} = currentToken;
                end
            end

            %% Nested function to parse a token into IDs
            function ids = parseToken(token)
                if startsWith(token, '[') && endsWith(token, ']')
                    % It's a group
                    content = token(2:end-1);
                    if isempty(content)
                        error('Empty group in connection string.');
                    end
                    idStrs = strsplit(content, ',');
                    try
                        ids = cellfun(@str2double, idStrs);
                    catch
                        error('Invalid IDs in connection string.');
                    end
                else
                    % It's a single ID
                    try
                        ids = str2double(token);
                    catch
                        error('Invalid ID in connection string.');
                    end
                end
            end
        end

        %% Add Stimuli to a Dendrite
        function addStimuli(obj, ID, Signal)
            if ~obj.isValidID(ID)
                error('Dendrite ID %d is invalid.', ID);
            end
            idx = obj.getDendriteIndex(ID);
            obj.Stimuli{idx} = Signal;
        end

        %% Verify Connectivity (All Dendrites Connected to Soma)
        function verify(obj)
            % Verify that all dendrites are connected to the soma (ID = 0)
            % If there are floating dendrites, output an error with their IDs

            % Get the list of all IDs
            allIDs = [obj.dendrites.ID];

            % Initialize visited set
            visitedIDs = [];

            % Initialize a queue for BFS
            queue = [];

            % Start from soma (ID = 0)
            if ~obj.isValidID(0)
                error('Soma (ID = 0) not found in the model.');
            end

            queue(end+1) = 0;
            visitedIDs(end+1) = 0;

            while ~isempty(queue)
                currentID = queue(1);
                queue(1) = []; % Remove the first element

                idx = obj.getDendriteIndex(currentID);
                conn = obj.Connectivity(idx);

                % Get connected IDs (both proximal and distal)
                connectedIDs = [conn.List_Proximal, conn.List_Distal];

                for i = 1:length(connectedIDs)
                    connectedID = connectedIDs(i);
                    if ~ismember(connectedID, visitedIDs)
                        visitedIDs(end+1) = connectedID;
                        queue(end+1) = connectedID;
                    end
                end
            end

            % Now compare visitedIDs with allIDs
            unvisitedIDs = setdiff(allIDs, visitedIDs);
            if ~isempty(unvisitedIDs)
                error('Floating dendrites detected with IDs: %s', mat2str(unvisitedIDs));
            end
            % If everything is connected, do nothing
        end
    end
end
