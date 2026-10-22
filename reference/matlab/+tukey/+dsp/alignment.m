function g = alignment()
    % ALIGNMENT Finding a template in a signal:
    % findsignal with its default options, as the target app calls it.
    g = struct("category", {"findsignal"}, "run", {@findsignalVector});
end

function [output, description] = findsignalVector(params)
    % findsignal gives 1-based istart and istop, written 0-based to match
    % SciPy. It gives no distance profile, only the best window's.
    [istart, istop, dist] = findsignal(params.data, params.signal);
    output = struct("start", istart - 1, "stop", istop - 1, ...
                    "dist", encode.vec(dist));
    description = "[istart, istop, dist] = findsignal(data, signal)";
end
