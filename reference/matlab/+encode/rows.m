function s = rows(m)
% rows encodes a matrix as a JSON array of its rows
% works on a single row
r = strings(size(m, 1), 1);
for i = 1:size(m, 1)
    r(i) = encode.vec(m(i, :));
end
s = "[" + strjoin(r, ", ") + "]";
end
