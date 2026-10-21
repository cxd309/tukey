function g = design()
    % DESIGN Filter design: butter.
    %
    % Not reproduced: zpk2sos, as zp2sos factorises differently (the same
    % filter in different sections).
    g = struct("category", {"butter"}, "run", {@butterVector});
end

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
