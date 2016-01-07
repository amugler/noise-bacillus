function D = JS_div(p, q)

D = -.5*(entropy(p) + entropy(q)) + entropy(.5*(p + q));
