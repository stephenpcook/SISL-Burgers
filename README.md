# Static and moving semi-Lagrangian solutions to Burgers' equation

## Code in Brief

This MATLAB code is designed for running simulations of Burgers' equation using
a semi-Lagrangian numerical scheme with a number of different interpolants and
time dependent meshes and estimating the front speed and front width.  It was
created by Stephen P. Cook between 2011 and 2018 to support research by Stephen P.
Cook, Chris Budd, Adrian Hill and Thomas Melvin.

Stephen P. Cook, 2018-01-03, s.cook@bath.ac.uk

## Basic Usage

These have been tested on windows with MATLAB 2016b.

### Initial Setup

Supporting files (.mat files) can be created with the function setup.m.

### Semi-Lagrangian Solution to Burgers' Equation

A single simulation of Burgers' equation can be run with the code in
burgersSLMM.m - see the file for usage notes.

### Numerical experiments with changing spatial and temporal steps

Sets of experiments were given alphabetic names: Alan and Barry. After the
initial setup, each individual experiment can be run with expt1.m - for
example, the first "Alan" experiment can be run with

    expt1('experiments/alan/alan1')

### Numerical experiments with fixed Courant numbers

An experiment option file can be run with fixed Courant numbers using the
function expt_cfl.m - for example, an experiment with a hermite interpolant can
be run with

    expt_cfl('options/interpolation/hermite')

## Maintenance

The author does not intend to maintain this code; it was not initially planned
to be distributed, but is made available to allow the results presented in the
research to be tractable.

## Layout of the code

The main functions are in burgersSLMM.m, expt1.m and expt_cfl.m, described
above under Basic Usage.  Also in the root folder are:
 - plotting functions (expt2.m, expt_cfl_plot.m and theory_plots.m)
 - functions used by burgersSLMM (feval_lim.m, fwd_euler.m and mmpde5.m)
 - a function for evaluating the front width and location (get_m_x.m)

Supporting .mat files created by setup.m are placed with their generating
functions in the subfolders in experiments and options.  These provide the
parameters to the solver.

Interpolation routines are found in the interpolation folder, and routines for
moving meshes are provided in the mm_suite folder.

The code is in various states of documentation and maintainability.

## License

This code is provided under the MIT license. See LICENSE.txt
