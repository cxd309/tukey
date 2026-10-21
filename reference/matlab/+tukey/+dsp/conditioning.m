function g = conditioning()
    % CONDITIONING Conditioning and reductions: detrending, interpolation,
    % filling gaps, RMS.
    categories = {"detrend", "interp1", "fillmissing", "rms"};
    runs = {@detrendVector, @interp1Vector, @fillmissingVector, @rmsVector};
    g = struct("category", categories, "run", runs);
end

function [output, description] = detrendVector(params)
    y = detrend(params.x, params.type);
    output = struct("y", encode.vec(y));
    description = sprintf("detrend(x, ""%s"")", params.type);
end

function [output, description] = interp1Vector(params)
    % linear by default, NaN outside the range of x
    yi = interp1(params.x, params.y, params.xi);
    output = struct("yi", encode.vec(yi));
    description = "interp1(x, y, xi)";
end

function [output, description] = fillmissingVector(params)
    if params.method == "constant"
        y = fillmissing(params.x, "constant", params.value);
        description = sprintf("fillmissing(x, ""constant"", %g)", params.value);
    else
        y = fillmissing(params.x, params.method);
        description = sprintf("fillmissing(x, ""%s"")", params.method);
    end
    output = struct("y", encode.vec(y));
end

function [output, description] = rmsVector(params)
    output = struct("rms", encode.vec(rms(params.x)));
    description = "rms(x)";
end
