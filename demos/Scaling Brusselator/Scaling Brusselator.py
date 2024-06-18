parameters={'sigma1':1,
            'sigma2':0.1,
            'delta':0.5,
            "xi": 1,
            "c": 4.8,
            "K1": 0.005,
            "K2": 1,
            "D1": 1,
            "D2": 1,
            "R": 1}

var=['u',
     'v']

diffmatrix=[['D1', 0, 0, 0],
            [0, 'D2', 0, 0],
            [0, 0, 'delta**2', 0],
            [0, 0, 0, 1]]

kinetics=["-sigma1*U",
        "-sigma2*V",
        "xi/delta - c*u + u**2*v - K1*(u - U)",
        "(c - 1)*u - u**2*v - K2*(v - V)"]

lvalues = [1, 2, 3, 4]

tol=1e-7

find_eq = 'n'

phiunit='n'

thirdcoef = 'y'

crosscoef = 'y'

crosspar = 'c'

plot2d='y'

parameters_on_axes = ['c', 'delta']

names_of_parameters = ['gamma', 'delta']

intervalx=[4, 15]

intervaly=[0, 1]

lines_to_search={'delta':0.5}