# radial-euler
[![Build Status](https://travis-ci.org/bvwells/radial-euler.svg?branch=master)](https://travis-ci.org/bvwells/radial-euler)

Compressible Euler Equations in radial coordinates solved with finite volume approach.

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
