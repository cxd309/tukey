function g = generators()
    % GENERATORS The tukey/dsp categories MATLAB can reproduce, and how
    % each run function takes a SciPy vector's params and returns MATLAB's
    % output (encoded with the encode namespace) and a description,
    % or an empty output when MATLAB has no equivalent for that case
    %
    % not reproduced:
    % lfilter_zi and sosfilt_zi (no MATLAB function)
    % zpk2sos (zp2sos factorises differently)
    categories = {"butter", "filter", "lfilter_state", "filtfilt", ...
                  "sosfilt", "sosfiltfilt", "freqz", "freqz_sos"};
    runs = {@butterVector, @filterVector, @filterStateVector, ...
            @filtfiltVector, @sosfiltVector, @sosfiltfiltVector, ...
            @freqzVector, @sosFreqzVector};
    g = struct("category", categories, "run", runs);
end

% -------------------------------------------------------------------------
% design
% -------------------------------------------------------------------------

function [output, description] = butterVector(params)
    % butter's (b, a), and its zpk through zp2sos for the sos output
    bandTypes = dictionary( ...
                           ["lowpass", "highpass", "bandpass", "bandstop"], ...
                           ["low", "high", "bandpass", "stop"]);
    btype = bandTypes(params.btype);
    wn = params.wn';

    [b, a] = butter(params.order, wn, btype);
    [z, p, k] = butter(params.order, wn, btype);
    sos = zp2sos(z, p, k);

    output = struct("b", encode.vec(b), "a", encode.vec(a), ...
                    "sos", encode.rows(sos));
    description = sprintf("butter(%d, %s, ""%s""); zp2sos of its zpk", ...
                          params.order, mat2str(wn), btype);
end

% -------------------------------------------------------------------------
% filtering
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% frequency response
% -------------------------------------------------------------------------

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
