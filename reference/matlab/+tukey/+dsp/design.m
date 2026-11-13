function g = design()
    % DESIGN Filter design: butter, cheby1, cheby2, ellip and bessel, the
    % order selections buttord, cheb1ord, cheb2ord and ellipord, fir1, and
    % iirnotch and iirpeak, which are in the DSP System Toolbox.
    %
    % Not reproduced: zpk2sos, as zp2sos factorises differently (the same
    % filter in different sections).
    categories = {"butter", "cheby1", "cheby2", "ellip", "bessel", ...
                  "buttord", "cheb1ord", "cheb2ord", "ellipord", "firwin", ...
                  "iirnotch", "iirpeak"};
    runs = {@butterVector, @cheby1Vector, @cheby2Vector, @ellipVector, ...
            @besselVector, ...
            @(params) ordVector(@buttord, "buttord", params), ...
            @(params) ordVector(@cheb1ord, "cheb1ord", params), ...
            @(params) ordVector(@cheb2ord, "cheb2ord", params), ...
            @(params) ordVector(@ellipord, "ellipord", params), ...
            @firwinVector, ...
            @(params) notchVector(@iirnotch, "iirnotch", params), ...
            @(params) notchVector(@iirpeak, "iirpeak", params)};
    g = struct("category", categories, "run", runs);
end

function [output, description] = ordVector(select, name, params)
    % An order selection's order and wn. MATLAB's buttord and cheb2ord meet
    % the stopband edge exactly, where SciPy meets the passband, and its
    % bandstop edges come from a numerical minimiser.
    wp = params.wp';
    ws = params.ws';
    [order, wn] = select(wp, ws, params.ripple_db, params.attenuation_db);
    output = struct("order", sprintf("%d", order), "wn", encode.vec(wn));
    description = sprintf("%s(%s, %s, %g, %g)", name, mat2str(wp), ...
                          mat2str(ws), params.ripple_db, params.attenuation_db);
end

function [output, description] = butterVector(params)
    [output, description] = designVector(@butter, "butter", params, {});
end

function [output, description] = cheby1Vector(params)
    [output, description] = designVector(@cheby1, "cheby1", params, ...
                                         {params.ripple_db});
end

function [output, description] = cheby2Vector(params)
    [output, description] = designVector(@cheby2, "cheby2", params, ...
                                         {params.attenuation_db});
end

function [output, description] = ellipVector(params)
    [output, description] = designVector(@ellip, "ellip", params, ...
                                         {params.ripple_db, ...
                                          params.attenuation_db});
end

function [output, description] = besselVector(params)
    % MATLAB has no digital Bessel design, so this is besself's analog design
    % at the prewarped edges, through bilinear at fs = 2, as SciPy's bessel.
    % besselap is normalised by phase, as SciPy's norm="phase", and stops at
    % order 25, so only those cases are reproduced.
    output = [];
    description = "";
    if params.norm ~= "phase" || params.order > 25
        return
    end
    bandTypes = dictionary( ...
                           ["lowpass", "highpass", "bandpass", "bandstop"], ...
                           ["low", "high", "bandpass", "stop"]);
    btype = bandTypes(params.btype);
    wn = params.wn';
    warped = 4 * tan(pi * wn / 2);

    [z, p, k] = besself(params.order, warped, btype);
    [z, p, k] = bilinear(z, p, k, 2);
    b = k * poly(z);
    a = poly(p);
    output = struct("b", encode.vec(b), "a", encode.vec(a), ...
                    "sos", encode.rows(zp2sos(z, p, k)));
    description = sprintf("besself(%d, 4*tan(pi*%s/2), ""%s""); " + ...
                          "bilinear at fs = 2", ...
                          params.order, mat2str(wn), btype);
end

function [output, description] = firwinVector(params)
    % fir1 of order numtaps - 1 with the case's window samples. fir1 needs
    % an order of at least 1, so the one-tap case is skipped.
    output = [];
    description = "";
    window = params.window';
    if numel(window) < 2
        return
    end
    bandTypes = dictionary( ...
                           ["lowpass", "highpass", "bandpass", "bandstop"], ...
                           ["low", "high", "bandpass", "stop"]);
    btype = bandTypes(params.btype);
    wn = params.wn';
    scaling = "scale";
    if ~params.scale
        scaling = "noscale";
    end

    h = fir1(numel(window) - 1, wn, btype, window, scaling);
    output = struct("h", encode.vec(h));
    description = sprintf("fir1(%d, %s, ""%s"", window, ""%s"")", ...
                          numel(window) - 1, mat2str(wn), btype, scaling);
end

function [output, description] = notchVector(design, name, params)
    % iirnotch or iirpeak, which take the -3 dB bandwidth w0/q, not q
    [b, a] = design(params.w0, params.w0 / params.q);
    output = struct("b", encode.vec(b), "a", encode.vec(a));
    description = sprintf("%s(%g, %g/%g)", name, params.w0, params.w0, ...
                          params.q);
end

function [output, description] = designVector(design, name, params, dbArgs)
    % A design's (b, a), and its zpk through zp2sos for the sos output.
    % dbArgs is the design's dB arguments, if it has any.
    bandTypes = dictionary( ...
                           ["lowpass", "highpass", "bandpass", "bandstop"], ...
                           ["low", "high", "bandpass", "stop"]);
    btype = bandTypes(params.btype);
    wn = params.wn';

    [b, a] = design(params.order, dbArgs{:}, wn, btype);
    [z, p, k] = design(params.order, dbArgs{:}, wn, btype);
    sos = zp2sos(z, p, k);

    output = struct("b", encode.vec(b), "a", encode.vec(a), ...
                    "sos", encode.rows(sos));
    dbText = "";
    if ~isempty(dbArgs)
        dbText = sprintf("%g, ", dbArgs{:});
    end
    description = sprintf("%s(%d, %s%s, ""%s""); zp2sos of its zpk", ...
                          name, params.order, dbText, mat2str(wn), btype);
end
