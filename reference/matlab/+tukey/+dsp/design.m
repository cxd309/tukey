function g = design()
    % DESIGN Filter design: butter, cheby1, cheby2 and ellip, and their
    % order selections buttord, cheb1ord, cheb2ord and ellipord.
    %
    % Not reproduced: zpk2sos, as zp2sos factorises differently (the same
    % filter in different sections).
    categories = {"butter", "cheby1", "cheby2", "ellip", ...
                  "buttord", "cheb1ord", "cheb2ord", "ellipord"};
    runs = {@butterVector, @cheby1Vector, @cheby2Vector, @ellipVector, ...
            @(params) ordVector(@buttord, "buttord", params), ...
            @(params) ordVector(@cheb1ord, "cheb1ord", params), ...
            @(params) ordVector(@cheb2ord, "cheb2ord", params), ...
            @(params) ordVector(@ellipord, "ellipord", params)};
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
