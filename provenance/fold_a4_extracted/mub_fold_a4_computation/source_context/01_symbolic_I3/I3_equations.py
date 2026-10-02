import sympy as sp
I=sp.I
z,t=sp.symbols("z t")
x=sp.symbols("x1:6"); y=sp.symbols("y1:6")
F=[
    x1*y1 + x1*y2 + x1*y3 + x1*y4 + x1*y5 + x1 + x2*y1 + x2*y2 + x2*y3 + x2*y4 + x2*y5 + x2 + x3*y1 + x3*y2 + x3*y3 + x3*y4 + x3*y5 + x3 + x4*y1 + x4*y2 + x4*y3 + x4*y4 + x4*y5 + x4 + x5*y1 + x5*y2 + x5*y3 + x5*y4 + x5*y5 + x5 + y1 + y2 + y3 + y4 + y5 - 5,
    -t*x1*y4 + t*x1*y5 + I*t*x2*y4 - I*t*x2*y5 - I*t*x3*y4 + I*t*x3*y5 + t*x4*y4*z - t*x4*y5*z - t*x5*y4*z + t*x5*y5*z + t*y4 - t*y5 + x1*y1 + I*x1*y2 - I*x1*y3 - x1 - I*x2*y1 + x2*y2 - x2*y3 + I*x2 + I*x3*y1 - x3*y2 + x3*y3 - I*x3 - x4*y1*z - I*x4*y2*z + I*x4*y3*z + x4*z + x5*y1*z + I*x5*y2*z - I*x5*y3*z - x5*z - y1 - I*y2 + I*y3 - 5,
    t*x1*y1*z + I*t*x1*y2 - t*x1*y3*z - I*t*x1*y4 - t*x1*y5 + t*x1 - t*x3*y1*z - I*t*x3*y2 + t*x3*y3*z + I*t*x3*y4 + t*x3*y5 - t*x3 - I*x2*y1*z + x2*y2 + I*x2*y3*z - x2*y4 + I*x2*y5 - I*x2 + I*x4*y1*z - x4*y2 - I*x4*y3*z + x4*y4 - I*x4*y5 + I*x4 - x5*y1*z - I*x5*y2 + x5*y3*z + I*x5*y4 + x5*y5 - x5 + y1*z + I*y2 - y3*z - I*y4 - y5 - 5,
    t*x1*y1*z - I*t*x1*y2 - t*x1*y3*z + t*x1*y4 + I*t*x1*y5 - t*x1 - t*x3*y1*z + I*t*x3*y2 + t*x3*y3*z - t*x3*y4 - I*t*x3*y5 + t*x3 + I*x2*y1*z + x2*y2 - I*x2*y3*z + I*x2*y4 - x2*y5 - I*x2 + x4*y1*z - I*x4*y2 - x4*y3*z + x4*y4 + I*x4*y5 - x4 - I*x5*y1*z - x5*y2 + I*x5*y3*z - I*x5*y4 + x5*y5 + I*x5 - y1*z + I*y2 + y3*z - y4 - I*y5 - 5,
    I*t*x1*y4 - I*t*x1*y5 - I*t*x2*y4 + I*t*x2*y5 + t*x3*y4 - t*x3*y5 + t*x4*y4*z - t*x4*y5*z - t*x5*y4*z + t*x5*y5*z - t*y4 + t*y5 + x1*y1 - x1*y2 + I*x1*y3 - I*x1 - x2*y1 + x2*y2 - I*x2*y3 + I*x2 - I*x3*y1 + I*x3*y2 + x3*y3 - x3 - I*x4*y1*z + I*x4*y2*z + x4*y3*z - x4*z + I*x5*y1*z - I*x5*y2*z - x5*y3*z + x5*z + I*y1 - I*y2 - y3 - 5,
    x1*y1 - I*x1*y2 + x1*y3 - x1*y4 - x1*y5 + I*x1 + I*x2*y1 + x2*y2 + I*x2*y3 - I*x2*y4 - I*x2*y5 - x2 + x3*y1 - I*x3*y2 + x3*y3 - x3*y4 - x3*y5 + I*x3 - x4*y1 + I*x4*y2 - x4*y3 + x4*y4 + x4*y5 - I*x4 - x5*y1 + I*x5*y2 - x5*y3 + x5*y4 + x5*y5 - I*x5 - I*y1 - y2 - I*y3 + I*y4 + I*y5 - 5,
]
