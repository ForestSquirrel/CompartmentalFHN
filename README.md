# Dendrite tree model simulator based on FirzHugh-Nagumo model

## Table of Contents
- [Dendrite tree model simulator based on FirzHugh-Nagumo model](#dendrite-tree-model-simulator-based-on-firzhugh-nagumo-model)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Features](#features)
  - [Project Structure](#project-structure)
  - [Examples](#examples)
  - [Experiments](#experiments)
- [Documentation](#documentation)
  - [Classes](#classes)
    - [Dendrite](#dendrite)
    - [DendriteTreeModel](#dendritetreemodel)
  - [Toolbox Functions](#toolbox-functions)
    - [treeSolver](#treesolver)
    - [classifyFiringPattern](#classifyfiringpattern)
  - [Utils Functions](#utils-functions)
    - [Chirp](#chirp)
    - [NFR](#nfr)
    - [supraStep](#suprastep)
  - [CPP Functions](#cpp-functions)
    - [mexRK4Solver](#mexrk4solver)


## Introduction
This project, the **Dendrite Tree Model Simulator based on the FitzHugh-Nagumo (FHN) model**, is a MATLAB-based framework designed to simulate and analyze the dynamics of dendritic tree structures using compartmental models.

The generalized equation for the FHN compartment, around which this project is structured, is as follows:

<img src="https://latex.codecogs.com/svg.image?\begin{cases}\frac{d&space;u}{d&space;t}=\tau(NaX\cdot(u(u-1)(1-\alpha&space;u)-v)&plus;\sum_{i=1}^{N_p}gc(u-u_{p_i})&plus;\sum_{i=1}^{N_d}gc(u-u_{d_i})&plus;I_{ext})\\\frac{d&space;v}{d&space;t}=\tau_r&space;b&space;u\end{cases}," title="FHN Compartment" />

where:
- $u$ - Membrane voltage
- $v$ - Gating current
- $NaX$ - Sodium channel density (e.g., $NaX = 1$ for soma)
- $\tau, \tau_r$ - Time constants
- $I_{ext}$ - External input
- $gc$ - Conductance between compartments
- $\sum_{i=1}^{N_p} gc(u - u_{p_i})$ and $\sum_{i=1}^{N_d} gc(u - u_{d_i})$ - Coupling terms from proximal $u_{p_i}$  and distal $u_{d_i}$  dendrites, with $N_p$  and $N_d$  as their respective counts
- Parameters $\alpha$  and $b$ - Control system nonlinearity

## Features
- **Modular Dendritic Tree Construction**: Build complex dendritic tree models by connecting dendrites to a central soma, allowing users to simulate the propagation of signals across various branch configurations.
- **Flexible Parameterization**: Users can configure properties like time constants, sodium channel density, and conductance, enabling a wide range of neuronal behavior modeling.
- **Dynamic Stimuli Application**: Apply various types of stimuli, such as step, sine, or chirp signals, to any part of the dendritic tree to observe the response patterns and study excitability.
- **Efficient Simulation with Runge-Kutta Integration**: A C++ MEX function enables fast numerical integration using the 4th-order Runge-Kutta method, allowing users to explore complex models with large numbers of dendrites efficiently.

## Project Structure
The project is stuctured as follows:
- [classes](/classes) - Folder with MATLAB classess
- [cpp](/cpp) - MEX stuff
  - [sources](/cpp/sources) - Source `.cpp` files
  - [compiled](/cpp/compiled) - Compiled `.mexw64`
- [experiments](/experiments) - MATALB live-scripts (`.mlx`) with experiments
- [toolbox](/toolbox) - MATLAB functions directly related to neural stuff
- [utils](/utils) - Utilitty MATALB functions, used in some experiments

> [!NOTE]  
> Also for those not familir with MATALB projects, this "internal" stuff is also included:
>
> - [CompatmentalFHN.prj](/CompatmentalFHN.prj) - Main project file, much like solution file (`.sln`) in generic coding.
> - [resources](/resources) - Internal MATLAB file mapings for project to work
>
> The main feature here is that all the `paths` will be automatically adjusted the moment you "run" `CompatmentalFHN.prj`. This way no dependencies will be absent from MATALB search path.

## Examples
- [example_V](/examples/example_V.m)
  
  Shows basic usage of the dendrite tree:
  1. Bulid a model with 2 parallel dendrites connected to soma
  2. Add sine stimuli to both of dendrites
  3. Plot soma response


## Experiments
- [BD](/experiments/BD.mlx) - Basic 1D bifurcation diagrams
- [FIR](/experiments/FIR.mlx) - Frequency response of a tree model
- [Misc](/experiments/Misc.mlx) - Testing of firing pattern classification based on [Hippocampone.org](https://hippocampome.org/) approach.

# Documentation
## Classes
### [Dendrite](/classes/Dendrite.m)
This class represents a single dendrite in the model. Each `Dendrite` instance has a unique ID, parameters for its behavior, and initial conditions for simulation.

**Properties:**
- `ID`: Unique identifier for the dendrite.
- `params`: Struct holding parameters:
  - `alpha`: Nonlinearity coefficient.
  - `b`: Nonlinearity coefficient.
  - `Tau`: Time constant for u.
  - `TauR`: Time constant for v.
  - `NaX`: Sodium channel density (e.g., 1 for soma).
  - `gc`: Conductance to other compartments.
- `InitialConditions`: Initial values for the model state `[u0; v0]`.

**Constructor:**
`Dendrite(ID, params, InitialConditions)`: Initializes a `Dendrite` object.
- `ID` (required): Unique identifier.
- `params` (optional): Parameters struct. If not provided, default values are used.
- `InitialConditions` (optional): Initial conditions `[u0; v0]`. If not provided, defaults to `[0.2; 0.2]`.
- 
### [DendriteTreeModel](/classes/DendriteTreeModel.m)
This class models a dendritic tree structure, managing an array of Dendrite objects, connections between them, and external stimuli. It supports operations to add dendrites, define connectivity, and ensure network integrity.

**Properties:**
- `dendrites`: Array of Dendrite objects in the model.
- `numDendrites`: Count of dendrites in the model.
- `Stimuli`: Cell array to hold stimuli associated with each dendrite.
- `Connectivity`: Struct array describing connections between dendrites, with fields for proximal and distal connections.
- `id_to_idx`: Mapping from dendrite ID to its index in dendrites.
- 
**Constructor:**

`DendriteTreeModel()`: Initializes an empty dendritic tree model.

**Methods:**
- `addDendrite(dendrite)`: Adds a Dendrite object to the model.
- `getDendriteIndex(ID)`: Retrieves the index of a dendrite by its ID.
- `isValidID(ID)`: Checks if a given dendrite ID is valid.
- `addConnection(ID, List_Proximal, List_Distal)`: Defines connections for a dendrite by specifying lists of proximal and distal connections.
- `addConnectionStr(connStr)`: Adds connections based on a formatted connection string.
- `addStimuli(ID, Signal)`: Assigns a stimulus to a specific dendrite.
- `verify()`: Validates that all dendrites are connected to the soma (ID = 0) to prevent isolated nodes.
## Toolbox Functions
### [treeSolver](/toolbox/treeSolver.m)
The `treeSolver` function acts as a wrapper, preparing input data from a `DendriteTreeModel` object and calling the `mexRK4Solver` MEX function for simulation. It returns the time vector and the solution matrix, containing the state variables for each dendrite at each time step.

**Inputs:**
- ***treeModel*** (`DendriteTreeModel` object): Model representing the dendritic tree.
- ***tmax*** (scalar): Maximum simulation time.
- ***h*** (scalar): Time step size.

**Outputs:**
- ***t*** (vector): Time vector.
- ***solution*** (matrix): Simulated state variables over time.

**Workflow:**
1. ***Prepare Time Vector***: Computes the time steps for the simulation.
2. ***Extract Parameters***: Collects dendrite parameters and connections from `treeModel`.
3. ***Initial Conditions***: Uses helper functions to initialize state variables.
4. ***Coupling and Stimulus Preparation***: Constructs connectivity data and stimulus matrix based on `treeModel`.
5. ***Run Simulation***: Calls `mexRK4Solver` with the prepared inputs and receives simulation results.

**Helper Functions:**
- `getInitialConditions`: Extracts initial conditions for each dendrite.
- `findCouplingMatrix`: Identifies connectivity and creates the coupling data matrix.
- `prepareStimuliMatrix`: Assembles stimuli data for the specified dendrites.
- `findStimulatedDendrites`: Finds dendrites with associated stimuli.
### [classifyFiringPattern](/toolbox/classifyFiringPattern.m)
The `classifyFiringPattern` function classifies the firing pattern of a neuron based on several metrics derived from a neuronal response signal, such as inter-spike intervals (ISI), first spike latency (FSL), post-spike silence (PSS), and the slow after hyperpolarizing wave amplitude (SWA). It follows the pseudocode methodology described in the [Hippocampome.org guide](https://hippocampome.org/php/Help_Firing_Pattern_Identification_Pseudocode.php).

**Inputs:**
- ***response*** (vector): The response signal of the neuron.
- ***fs*** (scalar): Sampling rate of the `response`.
- ***MinPeakHeight*** (scalar): Minimum peak height for identifying spikes.
- ***constants*** (optional, struct): Constants for the classification algorithm. Default values are used if not provided.

**Outputs:**
- ***firingPattern*** (string): A string describing the firing pattern, which can include markers for delay (`D.`), transient stuttering (`TSTUT.`), slow-wave bursting (`TSWB.`), steady state (`NASP/STEADY_STATE.`), persistent stuttering (`PSTUT.`), slow-wave bursting after spikes (`PSWB.`), and silence (`SLN.`).

**Steps and Classification Checks:**
1. Calculate ISI and First Spike Latency (FSL):
   - Identifies spike peaks and calculates ISI and FSL.
   - Checks if there are enough spikes to classify the pattern.
2. Compute Slow-Wave Amplitude (SWA) and Post-Spike Silence (PSS):
   - Uses the `computeSlowWave` subfunction to determine SWA and PSS based on identified spikes and inter-spike intervals.
3. Check for Delay:
   - Calls `hasDelay` to determine if the firing pattern includes a delayed response.
4. Check for Transient Stuttering or Slow-Wave Bursting:
   - Uses `hasTSTUT` to identify if there is transient stuttering or slow-wave bursting based on ISI patterns, SWA, and PSS.
5. Run Statistical Tests for Steady-State Firing:
   - Calls `runSolverStatTests` to fit various models and identify if the neuron reaches a steady state.
6. Check for Persistent Stuttering or Slow-Wave Bursting:
   - Calls `hasPSTUT` to identify persistent stuttering or slow-wave bursting.
7. Check for Silence:
   - Uses `hasSLN` to determine if there is silence after firing.

**Subfunctions:**
- `computeSlowWave:` Computes `swa` (slow-wave amplitude) and `pss` (post-spike silence) based on a slice of the response signal.
- `hasDelay`: Checks if the firing pattern exhibits delay.
- `hasTSTUT`: Identifies transient stuttering in the ISI pattern.
- `runSolverStatTests`: Performs piecewise linear fits to identify steady-state firing patterns.
- `hasPSTUT`: Checks for persistent stuttering in the firing pattern.
- `hasSLN`: Determines if the pattern has silence after firing.

## Utils Functions
### [Chirp](/utils/Chirp.m)
The `Chirp` function generates a frequency-modulated signal (chirp) with varying frequency over time, based on specified amplitude, frequency range, and chirp type.

**Inputs:**
- ***A*** (scalar): Amplitude of the signal. Must be a positive scalar.
- ***tmax*** (scalar): Duration of the signal in seconds. Must be a positive scalar.
- ***fs*** (scalar): Sampling rate in Hz. Must be a positive scalar.
- ***f0*** (scalar): Starting frequency in Hz. Must be a positive scalar.
- ***f1*** (scalar): Final frequency in Hz. Must be a positive scalar.
- ***options*** (struct, optional): Additional options for the chirp signal:
  - `type` (string, default: `'log'`): Specifies the type of chirp. Options:
    - `'linear'`: Linear frequency increase.
    - `'exp'`: Exponential frequency increase.
    - `'log'`: Logarithmic frequency increase.
  - ***signal*** (string, default: `'sin'`): Specifies the base waveform. Options:
    - `'sin'`: Sine waveform.
    - `'cos'`: Cosine waveform.
    - `'square'`: Square waveform.

**Outputs:**
- ***t*** (vector): Time vector for the duration of the signal.
- ***chirp_signal*** (vector): Generated chirp signal.
- ***f_t*** (vector): Instantaneous frequency at each time point.

**Chirp Types:**
- ***Linear***: Frequency varies linearly from `f0` to `f1` over `tmax`.
- ***Exponential***: Frequency increases exponentially from `f0` to `f1`.
- ***Logarithmic***: Frequency follows a logarithmic increase from `f0` to `f1`.

**Signal Types:**
- ***Sine***: Standard sine wave.
- ***Cosine***: Standard cosine wave.
- ***Square:*** Square wave.

**Example Usage:**
```matlab
[t, chirp_signal, f_t] = Chirp(1, 5, 1000, 20, 100, 'type', 'linear', 'signal', 'sin');
```
### [NFR](/utils/NFR.m)
The `NFR` (Nonlinear Frequency Response) function computes the frequency response of a system by comparing a chirp signal input and the corresponding output signal from the system. It leverages the Short-Time Fourier Transform (STFT) to estimate the system’s frequency response in terms of magnitude and phase across a range of frequencies over time.

**Inputs:**
- ***chirp_signal*** (vector): The input chirp signal.
- ***system_response*** (vector): The output signal of the system corresponding to the `chirp_signal `input.
- ***fs*** (scalar): Sampling rate of both chirp_signal and system_response.
- ***fft_opts*** (struct, optional): A struct for configuring FFT options with the following fields:
  - `window`: Window function for the STFT, default is `hamming(1024)`.
  - `nfft`: Number of FFT points. Default is the length of the `window`.
  - `noverlap`: Number of overlapping samples between segments. Default is 50% of `window` length.
  - `sides`: Specifies the type of FFT:
    - `'onesided'`: Default, returns the positive frequencies only.
    - `'twosided'`: Returns both positive and negative frequencies.
    - `'centered'`: Centers zero frequency.

**Outputs:**
- ***freq_response*** (struct): A struct containing the system's frequency response information:
- ***magnitude***: Magnitude response of the system.
- ***phase***: Phase response of the system.
- ***frequency***: Array of frequencies corresponding to the STFT.
- ***time***: Array of time instances corresponding to the STFT.
- 
**Example Usage:**
```matlab
% Define FFT options
fft_opts.window = hamming(512);
fft_opts.nfft = 512;
fft_opts.noverlap = 256;
fft_opts.sides = 'onesided';

% Compute the frequency response
freq_response = NFR(chirp_signal, system_response, fs, fft_opts);
```
**Helper Function:**

`getfield(struct, field, default)`: Retrieves the value of field from struct. If field does not exist, returns default.
### [supraStep](/utils/supraStep.m)
The `supraStep` function generates a customizable square wave signal with optional periods of zero amplitude at the beginning and end. It supports several types of waveforms and allows for signal inversion.

**Inputs:**
- ***Amplitude*** (scalar): Amplitude of the wave. Must be a positive scalar.
- ***Frequency*** (scalar): Frequency of the wave in Hz. Must be a positive scalar.
- ***SamplingRate*** (scalar): Sampling rate in Hz. Must be a positive scalar.
- ***NumPeriods*** (integer): Number of periods of the wave. Must be a positive integer.
- ***options*** (optional struct): A struct to configure additional waveform parameters:
  - `type` (string, default: `'square'`): Specifies the type of waveform:
    - `'square'`: Normal square wave ranging from `-Amplitude` to `Amplitude`.
    - `'step'`: Step wave ranging from `0` to `Amplitude`.
    - `'bipolar'`: Alternating waveform that cycles through `Amplitude`, `0`, `-Amplitude`, and `0`.
  - `inverse` (boolean, default: `false`): If true, inverts the signal.
  - `numZeroPeriodsStart` (integer, default: `1`): Number of zero-amplitude periods at the start.
  - `numZeroPeriodsEnd` (integer, default: `1`): Number of zero-amplitude periods at the end.

**Outputs:**
- ***t*** (vector): Time vector for the duration of the signal.
- ***sig*** (vector): Generated square wave signal.

**Waveform Types:**
- ***Square***: Standard square wave oscillating between `-Amplitude` and `Amplitude`.
- ***Step***: Step waveform oscillating between `0` and `Amplitude`.
- ***Bipolar***: Cycles through `Amplitude`, `0`, `-Amplitude`, and `0` within each period.

**Example Usage:**
```matlab
[t, sig] = supraStep(1, 1, 100, 2, 'type', 'bipolar', 'inverse', true, 'numZeroPeriodsStart', 2, 'numZeroPeriodsEnd', 1);
```

## CPP Functions
### [mexRK4Solver](/cpp/sources/mexRK4Solver.cpp)
This C++ MEX function performs the simulation of the dendritic tree model using the Runge-Kutta 4th-order (RK4) method. It calculates the evolution of state variables over time based on the parameters and connectivity of dendrites.

**Inputs:**

1. ***paramsArray*** (matrix, `numDendrites x 7`): Matrix of dendrite parameters.
    - Columns represent: `[alpha, b, Tau, TauR, NaX, gc, ID]`
2. ***couplingData*** (matrix, `numCouplings x 3`): Defines connectivity with `[row, col, value]`, where `row` and `col` are zero-based indices of connected dendrites, and `value` is the coupling conductance.
3. ***StimuliMatrix*** (matrix, `numStimDendrites x numTimeSteps`): Time series of stimuli applied to dendrites.
4. ***stimDendriteIDArray*** (vector): IDs of dendrites with applied stimuli, corresponding to rows of `StimuliMatrix`.
5. ***X0*** (vector): Initial conditions `[u0; v0]` for each dendrite.
6. ***h*** (scalar): Time step size.
7. ***numTimeSteps*** (scalar): Total number of time steps for simulation.
8. ***dendriteIDArray*** (vector): Array of dendrite IDs in the model.

**Outputs:**

1. ***t*** (vector): Time vector.
2. ***solution*** (matrix, `numStateVars x numTimeSteps`): Simulated state variables over time.
