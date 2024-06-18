# Amplitude-equations-bulk-surface-RD-equations
Scripts to compute the amplitude equations in a bulk-surface RD equation in a ball.

To get the coefficients of the amplitude equations, create a 'foo.py' file and place it into a folder called 'foo'. Said file must have the same format as the '.py' files in the 'demos' folder or 'name_of_folder.py'. Speficially, you need to provide the kinetics and the diffusion matrix of your system, together with names and values of parameters you would like to study in your bifurcation diagram.

After that, you run the Python file and the output will be some text files places into the folder 'foo', together with the Mathematica Plotter, which you can open and run to get the bifurcation curves and the values of the coefficients. Make sure you set the parameter values in said file as you would like to use them for the computation you are running.
