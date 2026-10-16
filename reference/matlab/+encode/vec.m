function s = vec(x)
% VEC Encode a real vector as a JSON array at full precision, 
% works on a single element
s = "[" + strjoin(compose("%.17g", x(:)'), ", ") + "]";
end
