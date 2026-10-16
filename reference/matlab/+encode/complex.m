function s = complex(z)
% COMPLEX Encode a complex vector as {"re": [...], "im": [...]}
s = "{""re"": " + encode.vec(real(z)) + ", ""im"": " + encode.vec(imag(z)) + "}";
end
