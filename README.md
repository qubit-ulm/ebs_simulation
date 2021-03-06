# Optical Tweezers Molecular Motor Simulation

This repository contains a simulations of molecular motor dynamics typically observed in optical tweezers. 
The simulation is written in matlab and contains simple poissonian stepping as well as
more complicated enzymatic dynamics like an elongating RNA Polymerase II (RNAP II). The latter example is motivated 
by recent experimental observations in Dangkulwanich et.al. eLife 2013;2:e00971

# Introduction

Dynamics of individual molecular motors can be directly observed by single molecule techniques such as optical tweezers.
Enzymatic movement happens in a step-like change of positions on a substrate. These steps can be as small as the distance of
two bases on DNA (0.34 nm). Since measurements have to be carried out at proper conditions for enzymes (room temperature and in liquids)  
steps are often highly corrupted by noise originating from brownian motion and the instrument itself. 
Therfore data simulation is a useful tool to assess further data analysis of these measurements. 
The following simulation methods can account for noise from harmonically trapped particles, instrument noise and can generate 
poissonian distributed steps as well as steps of an enzyme that switches between a translocation and paused state. The simulations use physical units. 
Positions are displayed in nm, Forces in pN and rates in Hz.

![Example Plot](https://github.com/qubit-ulm/pwcs_simulation/blob/master/examplestepimg.png)

# Usage

In matlab the following three commands generate noisy step data according to default parameters.

> sd = SimData();

> sd = sd.simulatePoissonSteps();

> sd = sd.simOptTrapNoise();

At first a SimData object is constructed. Then a class method that generates steps is called, here 'simulatePoissonSteps()'_ .
The array of step positions can be accessed through the 'SimData.pwcs'_ field, and e.g. plotted by:

> plot(sd.pwcs);

Finally noise is simulated around these steps using one of the noise simulation methods, here 'simOptTrapNoise()'_.
Noisy data can be accessed by the 'SimData.data'_ fields. There is also a 'plotData()'_ method that opens a figure and plots noisy 
as well as pure step data. Simply call:

> sd.plotData();

## General Remarks

A SimParams object contains all simulation parameters and the SimData class 
contains the simulated data arrays, step simulation methods and noise 
simulation methods. Noisy step data is simulated in three steps. First determine
all necessary parameters (like e.g. the number of data points) by creating a SimParams object. 
Note that there are already default parameters set when a SimParams object is created. 
Since these parameters are public members they can be directly accessed and changed. 
Secondly, a SimData object is created and noisy step data is generated by subsequent call of
a step generation method and finally a noise generation method.
The internals folder has to be added to matlab path.

# Further Example

The following commands simulate the movement of a Polymerase on a DNA of fixed length L pulled by a constant Force.

At first we set the simulation parameters:
> sp = SimParams();

> sp.N = 1e5;

> sp.h = 5e-4;

> sp.p2Pars.c_NTP = 500;

These commands produce similar noisy data as displayed in the figure above.
> sd = SimData(sp);

> sd = sd.simulatePol2();

> sd = sd.simOptTweezersVariableNoise();

We set the number of data points in sp.N to 1e5, and the time increments of the simulation in sp.h to 5e-4s.
'SimParams'_ contains an object 'p2Pars'_ that encapsulates all Polymerase relevant parameters. 
Here we have set an NTP concentration in p2Pars.c_NTP to 500 mM. In contrast to the upper example here the Polymerase
type simulation 'SimData.simulatePol2()'_ is used to generate steps. Moreover the noise simulation method
'SimData.simOptTweezersVariableNoise()'_ accounts for brownian motion in optical traps, instrumental fluctuations and 
the contribution of a successively shorter DNA tether.


# Documentation

## SimData class

### public attributes:

time

>time axis array

data

>noisy data array

pwcs

>a priori simulated steps array (without noise)

p

>time points of jumps, necessary for generating noisy data 

Simparams

>simulation parameters, SimParams object

### private attributes:

x_ins

>array of instrumental noise which was added to data

SimulationMethodSteps

>String of step generation method previously applied

SimulationMethodNoise

>String of noise generation method previously applied

### public member functions:

SimData

>constructor: 'sd = SimParams()'_ with default parameters or 'sd = SimParams(sp)'_                      

simulatePoissonSteps 

>simulates a composition of poisson process steps by using ./internals/simMultiplePoisson.m . Call
'sd = sd.simulatePoissonSteps()'_ .
Stepping rate is determined by 'SimData.Simparams.poissonPars'_ parameters. The number of poisson processes
is determined by the number of components of 'PoissonPars.jump_rate'_ , 'PoissonPars.delta'_ and 'PoissonPars.sigma_delta'_ .

simulatePol2

>simulates a transcribing Polymerase that can switch between elongation and paused state under constant force (specified by 'sd.Simparams.p2Pars.F'_ ) . Necessary parameters are stored 
in 'SimData.Simparams.p2Pars'_ . Call 'sd = sd.simulatePol2()'_ or 'sd = sd.simulatePol2(k1,kb1,kf,kb)'_ if specific elongation rate k1, 
entry to pause rate kb1, forward diffusion kf and backward diffusion kb are used. ./internals/simSimplePol2steps.m is used.

simulatePol2variableForce 

>simulates a transcribing Polymerase that can switch between elongation and paused state under changing external force. 
Initial force is specified by 'sd.Simparams.p2Pars.F'_ . 'F<0'_ corresponds to an opposing force, 'F>0'_ to an assisting force.
During stepping forces and transition rates are updated. 
Call 'sd = sd.simulatePol2variableForce()'_ . ./internals/simSimplePol2stepsincForce.m is used.

simOptTrapNoise   

>simulates noise which follows an Ornstein-Uhlenbeck process around the steps. Uses ./internals/simNoisyData.m and the optical trap simulation
parameters in 'sd.Simparams'_ . Call 'sd = sd.simOptTrapNoise()'_ , 'sd = sd.simOptTrapNoise(noiseFactor)'_ if vaiance of gaussian white noise should be increased by noiseFactor

simOptTweezersVariableNoise

>simulate noise according to a tethered bead system (two beads in optical traps tethered by DNA) with changing tether length. 
'sd = sd.simOptTweezersVariableNoise()'_ and 'sd = sd.simOptTweezersVariableNoise(noisefactor)'_ if vaiance of gaussian white noise should be increased by noiseFactor

plotData 

>plots data. 'sd.plotData()'_ 

getInstrNoise 

>returns an array of instrumental noise. Call 'y = sd.getInstrNoise()'_              

reScaleJumpInput 

>adapts data size when sd.Simparams.h and sd.Simparams.N are changed. '[pwcsneu, pneu] = sd.reScaleJumpInput(pwcs_old, p_old)'_              

whichSimulation

>returns a struc o two String describing . 'answ = sd.whichSimulation()'_


## SimParams class

### public attributes:

N 

>number of data points

h 

> time increments in [s]

kbT 

>thermal energy, room temperature in [pN*nm]

gamma 

> in viscosity of trapped bead in [pN*s/nm]

kappa 

>trap stiffness for trap1 and 2 [pN/nm]

k_DNA

>stiffness of a 1.8kbp long DNA tether (L=600nm)

L

>contour length of DNA in [nm]

x0

>initial elongtion from trap

k

>effective trap stiffness (kappa +2*k_DNA) of relative coordinate x

instrDrift

> correlation (0.1s) of instrument fluctuations in [1/s]

instrNoise

> standard deviation of fluctuations in [nm^2/s]

mix

>factor how much the beads signal is distorted by instr. Noise

simtype

> assist or oppose

constForce

> constant Force = true

p2Pars

> Pol2Params object

poissonPars

> PoissonParams object 

## PoissonParams class

public attributes:

jump_rate

> array of rates of poisson steps in [1/s]

delta

>array of stepsizes [nm]

sigma_delta

>standard deviation of step size in [nm]

## Pol2Params class

### public attributes

kbT

>thermal energy to calculate Boltzmann factors (e.g.: exp(F*delta/kbT) )of weights k_1, k_rev1, k_b1, k_bn, k_f

delta

>distance to transition state (0.17nm)

F

>external force in [pN]

k0

>unbiasd forward backward diffusion rate (at F=0)

k_f

>forward diffusion rate biased due to external force F

k_bn

>backward diffusion rate biased due to external force F

k_b1

>pause entry rate biased due to external force F

kb0

>unbiased pause entry rate (at F=0)

K_D

>NTP insertion equilibrium constant 

c_NTP

> NTP concentration in [mM]

k_1

> forward translocation

k_rev1

> reverse translocation

k_10

> forward translocation (F=0)

k_rev10

>reverse translocation (F=0)

k_3

>NTP catalysis rate

k2_net

> net NTP binding-catalysis rate. Approximates NTP binding + catalsis as a single step

k1_net

> net forward translocation-NTP binding-catalysis rate

k1

> net forward elongation (1/k1_net + 1/k2_net + 1/k_3)^-1


# License
The pwcs project is licensed to you under the Apache License, Version 2.0 (the "License"); you may not use the code except in compliance with the License. You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
