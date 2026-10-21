function s = vec(x)
    % VEC Encode a real vector as a JSON array at full precision, even with
    % one element. NaN becomes null, as JSON has no NaN
    s = "[" + strjoin(compose("%.17g", x(:)'), ", ") + "]";
    s = replace(s, "NaN", "null");
end
