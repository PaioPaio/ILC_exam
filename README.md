# ILC_exam

This repository attempts to reproduce the results of the paper "Iterative Learning in Functional Space for Non-Square Linear Systems" of C. Della Santina and F. Angelini. ([IEEE link here](https://ieeexplore.ieee.org/document/9683673)).


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia (install either from your distribution packet manager or [here](https://julialang.org/downloads/)) console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install and compile all necessary packages (might take a while) for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

TIP: You may want to start Julia with more than one thread to speed up the compilation process (and the code runs faster :O). To do that you can add the `-t auto` flag to your `julia` command (i.e. `julia -t auto` to start julia in the terminal).