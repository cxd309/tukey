function g = filtering()
    % FILTERING Filtering with (b, a) and second-order sections: at rest,
    % from a state, and zero-phase.
    %
    % Not reproduced: lfilter_zi and sosfilt_zi (no MATLAB function).
    categories = {"filter", "lfilter_state", "filtfilt", ...
                  "sosfilt", "sosfiltfilt", "sosfiltfilt_by_section", ...
                  "filtfilt_matlab"};
    runs = {@filterVector, @filterStateVector, @filtfiltVector, ...
            @sosfiltVector, @sosfiltfiltVector, @sosfiltfiltBySectionVector, ...
            @filtfiltMatlabVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = filterVector(params)
    y = filter(params.b, params.a, params.x);
    output = struct("y", encode.vec(y));
    description = "filter(b, a, x)";
end

function [output, description] = filterStateVector(params)
    [y, zf] = filter(params.b, params.a, params.x, params.zi);
    output = struct("y", encode.vec(y), "zf", encode.vec(zf));
    description = "[y, zf] = filter(b, a, x, zi)";
end

function [output, description] = filtfiltVector(params)
    % MATLAB has one edge handling, odd extension by 3 * filter order,
    % so only SciPy's default case has an equivalent
    output = [];
    description = "";
    if ~isDefaultPadding(params)
        return
    end
    y = filtfilt(params.b, params.a, params.x);
    output = struct("y", encode.vec(y));
    description = "filtfilt(b, a, x)";
end

function [output, description] = sosfiltVector(params)
    % sosfilt takes no initial state, so the with-state cases have no equivalent
    output = [];
    description = "";
    if ~isempty(params.zi)
        return
    end
    y = sosfilt(params.sos, params.x);
    output = struct("y", encode.vec(y), "zf", "null");
    description = "sosfilt(sos, x)";
end

function [output, description] = sosfiltfiltVector(params)
    % as filtfiltVector: only SciPy's default edge handling has an equivalent
    output = [];
    description = "";
    if ~isDefaultPadding(params)
        return
    end
    y = filtfilt(params.sos, 1, params.x);
    output = struct("y", encode.vec(y));
    description = "filtfilt(sos, 1, x)";
end

function [output, description] = sosfiltfiltBySectionVector(params)
    % filtfilt on sections filters each in turn, padded by 3 * the total
    % order, so only the cases with MATLAB's padding have an equivalent.
    % The gain is in the first section already, so g is 1.
    output = [];
    description = "";
    if ~params.matlab_padlen
        return
    end
    % a single section is a 1-by-6 row, which filtfilt warns could be b;
    % its fourth coefficient, 1, makes it read it as a section, as meant
    state = warning("off", "signal:filtfilt:ParseSOS");
    restore = onCleanup(@() warning(state));
    y = filtfilt(params.sos, 1, params.x);
    output = struct("y", encode.vec(y));
    description = "filtfilt(sos, 1, x)";
end

function [output, description] = filtfiltMatlabVector(params)
    % (b, a) with trailing zeros, whose order filtfilt counts without them
    y = filtfilt(params.b, params.a, params.x);
    output = struct("y", encode.vec(y));
    description = "filtfilt(b, a, x)";
end

function isDefault = isDefaultPadding(params)
    % SciPy's default edge handling: odd extension, padlen left to the default
    isDefault = isequal(params.padtype, 'odd') && isempty(params.padlen);
end
