# ILC_exam

This repository reproduces the results of the paper "Iterative Learning in Functional Space for Non-Square Linear Systems" of C. Della Santina and F. Angelini. ([IEEE link here](https://ieeexplore.ieee.org/document/9683673)).


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named

To (locally) reproduce this project, do the following:

0. Clone this code base.
1. Install Julia, either from your distribution packet manager or [here](https://julialang.org/downloads/) (and add it to your PATH if you're on Windows).
2. Open a terminal in the cloned folder and type "julia" to start a Julia REPL
3. Type "]" to enter package mode and download and compile the project dependencies as such:
   ```
   julia> ] activate .
   julia> ] instantiate
   ```

This will install and compile all necessary packages for you to be able to run the scripts and everything should work out of the box. While Julia is fast to install, this might take a while as you're practically installing the Julia version of Simulink.

*TIP:* You may want to start Julia with more than one thread to speed up the compilation process (and the code runs faster :O). To do that you can add the `-t auto` flag to your `julia` command (i.e. `julia -t auto` to start julia in the terminal).

4. To then play with the script in a MATLAB-like experience, use [VSCode](https://code.visualstudio.com/) and install the Julia extension. You can change the aforementioned number of threads in the extension setting to automatically launch Julia with a set number of threads.
5. Open the project inside the IDE and open either "massspringdamper.jl" or "basketball.jl" in the src folder.
6. Send the code blocks (delimited by "##") to the interactive terminal (REPL) with ALT/CTRL+Enter.
