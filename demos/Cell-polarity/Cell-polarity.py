parameters={"sigma1": 100.0,
            "sigma2": 1.0,
            "sigma3": 100.0,
            "sigma4": 1.0,
            "D1": 1,
            "D2": 1,
            "D3": 1,
            "D4": 1,
            "alpha": 0.05,
            "epsilon": 0.76,
            "mu": 0.5,
            "rho": 0.06,
            "xi": 2.7,
            "theta": 5.5,
            "eta": 3.6,
            "K1": 1.0,
            "K2": 100.0,
            "K3": 1.0,
            "K4": 100.0,
            "delta1": 0.33,
            "delta2": 1.0,
            "R": 1.0}

var=['u',
     'v',
     'w',
     'z']

diffmatrix=[["D1", "0", "0", "0", "0", "0", "0", "0"],
            ["0", "D2", "0", "0", "0", "0", "0", "0"],
            ["0", "0", "D3", "0", "0", "0", "0", "0"],
            ["0", "0", "0", "D4", "0", "0", "0", "0"],
            ["0", "0", "0", "0", "delta1**2", "0", "0", "0"],
            ["0", "0", "0", "0", "0", "1", "0", "0"],
            ["0", "0", "0", "0", "0", "0", "delta2**2", "0"],
            ["0", "0", "0", "0", "0", "0", "0", "1"]]

kinetics=["-sigma1*U",
        "-sigma2*V",
        "-sigma3*W",
        "-sigma4*Z",
        "(eta*u**2 + rho)*v - (alpha*w + mu)*u - epsilon*xi*u - K1*(u - U)",
        "-((eta*u**2 + rho)*v - (alpha*w + mu)*u) + epsilon*theta - K2*(v - V)",
        "(eta*w**2 + rho)*z - (alpha*u + mu)*w - epsilon*xi*w - K3*(w - W)",
        "-((eta*w**2 + rho)*z - (alpha*u + mu)*w) + epsilon*theta - K4*(z - Z)"]

lvalues = [1, 2, 3, 4]

tol = 1e-7

find_eq = 'n'

phiunit = 'n'

thirdcoef = 'y'

crosscoef = 'y'

crosspar = 'epsilon'

equilibrium = []

plot2d = 'y'

parameters_on_axes = ['epsilon', 'delta1']

names_of_parameters = ['kappa', 'delta1']

intervalx = [0.0, 1.0]

intervaly = [0, 0.8]

lines_to_search = {'delta1': 0.02, 'epsilon': 0.6}