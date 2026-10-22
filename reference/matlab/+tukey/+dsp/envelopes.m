function g = envelopes()
    % ENVELOPES The analytic signal and envelopes: hilbert, and envelope in
    % its analytic and RMS modes. Named envelopes, not envelope, so a call
    % to MATLAB's envelope in here doesn't read as this file calling itself.
    g = struct("category", {"hilbert", "envelope"}, ...
               "run", {@hilbertVector, @envelopeVector});
end

function [output, description] = hilbertVector(params)
    output = struct("analytic", encode.complex(hilbert(params.x)));
    description = "hilbert(x)";
end

function [output, description] = envelopeVector(params)
    if params.mode == "rms"
        [up, lo] = envelope(params.x, params.n, "rms");
        description = sprintf("[up, lo] = envelope(x, %d, ""rms"")", ...
                              params.n);
    else
        [up, lo] = envelope(params.x);
        description = "[up, lo] = envelope(x)";
    end
    output = struct("upper", encode.vec(up), "lower", encode.vec(lo));
end
