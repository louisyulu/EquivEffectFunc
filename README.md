# README #

This package implements the equivalent effect function (EEF) and intransic mode function (IMF) described at [Equivalent Effect Function and Fast Intrinsic Mode Decomposition](https://www.semanticscholar.org/paper/Equivalent-Effect-Function-and-Fast-Intrinsic-Mode-Lu/b8e35a90ce840041c4101ce28a19e6e223ab7898) by the author.

The algorithm uses polynomial spline function to aproximate the one dimensional data (or time series data). It is implemented in Julia language and depends on Dierckx package to handle the polynomial spline from order 1 to 5. It uses Pluto and Plotly packages to display the example results.


Here is a brief description of each implentation file:
- `EquivEffectFunc.jl` - EEF with control points depending on extrema and on change rate densities
- `control_points.jl` - two method of finding control points: 1) on consecutive maxima and minima of the original data value, first derivative or second derivative; 2) partition on absolute value (change rate) of original data value, first derivative or second derivative.
- `calculus.jl` - methods to calculate the derivatives and integral discretely.
- `examples_nb.jl` - demonstrates randomly generated data and real historical stock data displayed with Pluto notebook and Plotly charting library.

### Setup steps: ###
1. `git clone https://github.com/louisyulu/EquivEffectFunc.git EquivEffectFunc` 
Clone the repository
2. `cd EquivEffectFunc` 
Change to the package folder
3. `julia --project=.`
Start julia and activate the current project
4. `]` 
Type ] key to switch into package mode
5. `build`
Build the current package
6. Press `delete` key to return to Julia mode
7. `import Pluto`
Import Pluto package
8. `Pluto.run()`
Start Pluto notebook and the notebook will open in your default browser
9. Choose `examples_nb.jl` to start your experimentation. Have fun!
