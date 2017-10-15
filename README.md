# radial-euler
[![Build Status](https://travis-ci.org/bvwells/radial-euler.svg?branch=master)](https://travis-ci.org/bvwells/radial-euler)

Compressible Euler Equations in radial coordinates solved with finite volume approach.

The compressible Euler equations in radial coordinates are described by the non-linear system of partial differential equations

```
ρ_t + (ρu)_r = -(d-1)ρu/r
(ρu)_t + (ρu^2 + P)_r = -(d-1)ρu^2/r
E_t + (u(E+P))_r = -(d-1)u(E+P)/r
```

where ```ρ```, ```u```, ```P``` and ```E``` are the density, velocity, pressure and energy of the gas being modelled. The parameter ```d``` denotes the dimension of the problem. Values from one to three equate to cartesian, cylindrical and sphericl coordinates respectively. The equation of state for an ideal gas is used to close the system and is given by
```
E=P/(1-γ) + 0.5ρv^2
```
where ```γ``` is the ratio of specific heats for the gas.

## Building and Developing

Developing locally requires Docker for Windows. Run the command

```
docker build -t radial-euler .
```

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t -v /f/git/src/github.com/bvwells/radial-euler:/app radial-euler
```

This command maps the local workspace into the running image so any changes made in the running image will be reflected on the local workspace.

Within the running image generate the make files for the release version by running the command:

```
cmake .
```

To build the debug version of the code run the command:

```
cmake -DCMAKE_BUILD_TYPE=Debug
```

Build the executable by running the command:

```
make
```

## Running

The program takes the file [variables.data](./variables.data) as input to the simulation. The program can be run from the base of the repo with the command:

```
./bin/radialeuler.exe
```

The program outputs the solution into the files ```exact.m```.

## Plotting Solution

The output from the simulation can be plotted in [Octave](https://www.gnu.org/software/octave/) by running the plotting file
[plot_solution.m](./plot_solution.m) in the root of the repo.
