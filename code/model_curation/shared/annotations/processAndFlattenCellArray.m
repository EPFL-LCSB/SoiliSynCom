function flattenedCellArray = processAndFlattenCellArray(originalCellArray)
    % Initialize new cell array with the same size
    newCellArray = cell(size(originalCellArray));

    for i = 1:length(originalCellArray)
        % Check if the cell contains 'Not found'
        if ischar(originalCellArray{i}) || isstring(originalCellArray{i})
            if strcmp(originalCellArray{i}, 'Not found')
                newCellArray{i,1} = ''; % Assign a visibly empty cell
            else
                newCellArray{i,1} = originalCellArray{i}; % Copy as is
            end
        elseif iscell(originalCellArray{i})
            % Specific case for cell arrays with two or more elements
            if numel(originalCellArray{i}) >= 2
                newCellArray{i, 1} = originalCellArray{i}(end); % Keep only the last element
            elseif isempty(originalCellArray{i})
                newCellArray{i, 1} = ''; % Keep the cell visibly empty if it's already empty
            else
                % For cell arrays with one element or more (excluding 2x1 case)
                tempArray = originalCellArray{i};
                for j = 1:length(tempArray)
                    if strcmp(tempArray{j}, 'Not found')
                        tempArray{j} = ''; % Replace 'Not found' with a visibly empty cell
                    end
                end
                newCellArray{i, 1} = tempArray;
            end
        else
            % For any other type or condition not explicitly mentioned
            newCellArray{i, 1} = originalCellArray{i};
        end
    end

    % Initialize a new cell array for flattened content
    flattenedCellArray = cell(size(newCellArray));

    for i = 1:length(newCellArray)
        if iscell(newCellArray{i})
            % Check if it's a nested cell and not empty
            if ~isempty(newCellArray{i})
                % Extract the first element if it's nested
                flattenedCellArray{i} = newCellArray{i}{1};
            else
                % Handle empty nested cells or keep as is based on your requirement
                flattenedCellArray{i} = [];
            end
        else
            % Directly copy non-cell elements (handles 0x0 double as is)
            flattenedCellArray{i} = newCellArray{i};
        end
    end
end
