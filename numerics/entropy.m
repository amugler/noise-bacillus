function H = entropy(p)

p = p(find(p)); % removes zeros
H = -sum(p.*log2(p)); % in bits