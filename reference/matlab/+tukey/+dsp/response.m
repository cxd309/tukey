function g = response()
    % RESPONSE Frequency response and group delay of (b, a) filters and
    % second-order sections. Named response, not freqz, so a call to
    % MATLAB's freqz in here doesn't read as this file calling itself.
    g = struct("category", {"freqz", "freqz_sos", ...
                            "group_delay", "group_delay_sos"}, ...
               "run", {@freqzVector, @sosFreqzVector, ...
                       @groupDelayVector, @sosGroupDelayVector});
end

function [output, description] = freqzVector(params)
    % freqz needs at least one point,
    % so SciPy's empty n = 0 case has no equivalent
    output = [];
    description = "";
    if params.n == 0
        return
    end
    [h, w] = freqz(params.b, params.a, params.n);
    output = struct("w", encode.vec(w), "h", encode.complex(h));
    description = sprintf("[h, w] = freqz(b, a, %d)", params.n);
end

function [output, description] = sosFreqzVector(params)
    [h, w] = freqz(params.sos, params.n);
    output = struct("w", encode.vec(w), "h", encode.complex(h));
    description = sprintf("[h, w] = freqz(sos, %d)", params.n);
end

function [output, description] = groupDelayVector(params)
    % grpdelay needs at least one point, as freqz does
    output = [];
    description = "";
    if params.n == 0
        return
    end
    [delay, w] = grpdelay(params.b, params.a, params.n);
    output = struct("w", encode.vec(w), "delay", encode.vec(delay));
    description = sprintf("[gd, w] = grpdelay(b, a, %d)", params.n);
end

function [output, description] = sosGroupDelayVector(params)
    [delay, w] = grpdelay(params.sos, params.n);
    output = struct("w", encode.vec(w), "delay", encode.vec(delay));
    description = sprintf("[gd, w] = grpdelay(sos, %d)", params.n);
end
