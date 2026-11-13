function g = smoothing()
    % SMOOTHING Smoothing: the median filter.
    categories = {"medfilt"};
    runs = {@medfiltVector};
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
