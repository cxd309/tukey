function g = smoothing()
    % SMOOTHING Smoothing: the median filter and the Savitzky-Golay filter.
    categories = {"medfilt", "savgol_filter"};
    runs = {@medfiltVector, @savgolVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = medfiltVector(params)
    % 'truncate' takes the median of the samples inside x, and 'zeropad'
    % (the default) pads with zeros. A window holding NaN gives NaN, as the
    % default 'includenan'.
    y = medfilt1(params.x, params.n, char(params.edge));
    output = struct("y", encode.vec(y));
    description = sprintf("medfilt1(x, %d, ""%s"")", params.n, params.edge);
end

function [output, description] = savgolVector(params)
    % sgolayfilt fits the first and last framelen samples for the ends, as
    % SciPy's mode='interp'. It only smooths: a derivative has no equivalent.
    if params.deriv ~= 0
        output = [];
        description = "";
        return
    end
    y = sgolayfilt(params.x, params.polyorder, params.n);
    output = struct("y", encode.vec(y));
    description = sprintf("sgolayfilt(x, %d, %d)", params.polyorder, params.n);
end
