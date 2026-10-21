function g = filtering()
    % FILTERING Filtering with (b, a) and second-order sections: at rest,
    % from a state, and zero-phase.
    %
    % Not reproduced: lfilter_zi and sosfilt_zi (no MATLAB function).
    categories = {"filter", "lfilter_state", "filtfilt", ...
                  "sosfilt", "sosfiltfilt"};
    runs = {@filterVector, @filterStateVector, @filtfiltVector, ...
            @sosfiltVector, @sosfiltfiltVector};
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

function isDefault = isDefaultPadding(params)
    % SciPy's default edge handling: odd extension, padlen left to the default
    isDefault = isequal(params.padtype, 'odd') && isempty(params.padlen);
end
