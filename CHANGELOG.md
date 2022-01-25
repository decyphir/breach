# Release 1.10.0
- updated README.md
- **really** updated README.md
- like, now, it might be worth reading
- like, it even has an example of BreachSystem with a custom sim function
- and more details about BreachTraceSystem for monitoring data
- though there had to be some fixes in those two classes
- also minor changes in demo files
- plus minor things I don't remember

# Release 1.9.1

- Fixed demo BrDemo.from_file_signal_gen
- BreachTraceSystem new method GenTraceIdx

# Release 1.9.0 

- New solver: snobfit
- New additive semantics, to be documented
- Added Gui framework BreachGuiClass
- various bugfixes and performance improvements 

# Release 1.8.0

- Added past operators once and historically
- Better handling of corners in falsification
- Added binaries for macOS
- Various bugfixes

# Release 1.7.0

- binaries added removes need to run InstallBreach for windows and linux
- added ComputeMorrisSensi for sensitivity analysis
- Support for IA-STL 
- Various improvement in BreachProblem
- Simple pseudo random solver 

# Release 1.6.2

- Enveloppe computation in signal_gen_gui works with constraints
- New Global Nelder Mead properly re-established as default solver 
- Bugfix: precondition with non-simulink systems working again

# Release 1.6.1

- set domain cfg only does warning when trying to set a non existent parameter 

# Release 1.6.0

- Performance and memory overhead improvement when logging during falsification
- Global nelder mead default solver has been revamped 
- Diagnosis and IO-aware STL support added (see HSCC'19 paper)
 

## Bugfixes 
- Signature computation sometimes failed to properly link aliases 

# Release 1.5.2

- minor fix for maxstep/minstep detection causing crash with some solver configuration (Simulink)

# Release 1.5.1

- better detection of maximum and minimum step sizes when interfacing
  a model 
- some cleaning in Examples folder 

## Bugfixes
- fixed bug in SetInputGenGui 
- fixed bug whereas time was sometimes not set properly with SetTime
- fixed crash when minimum step size is bigger than eps, used for
  detecting logged signals when creating interface
- fixed bug causing BreachSimulinkSystem to fail when a logged object
  is not a simple signal 
- fixed bug when sampling 1D grids using legacy Refine 
- fixed typo bug in InstallBreach

# Release 1.5.0

- from_file_signal_gen can have an initialization script 
- mixed integer optimization support
- signal_gen_gui more robust 
- param_gen framework: arbitrary parameter transformation 
- better backward compatibility between BreachRequirement and STL_Formula 
- BreachSamplesPlot and BreachSignalsPlot improved
- BreachSimulinkSystem now look at all scalar parameters in the base
  workspace for potential tunable parameters, regardless if they are
  detected in the model
- STL formulas with no signals now are accepted (constraints on parameters)

## Bugfixes
- spike and exponential signal generators fixed
- fixed an issue with aliases not being always found 
- fixed duplicates in signals map

# Release 1.4.3

## Bugfixes
- display result fixed when using parallel falsification
- parallel computation during Nelder Mead phase disabled until
  behavior is fixed
  
# Release 1.4.2 

- Improved display values for multi-objective optimization

# Release 1.4.1

## Users 
- Simple support for multi-objective falsification with global Nelder-Mead 
- Added small example with Automatic transmission

## Bugfix
- traj_ref support fixed for BreachRequirement, would cause to fetch wrong 
  signal values for requirement evaluation some times  
- small bug in BreachSamplesPlot when instantiating contextual menus 

# Release 1.4.0

## Users
- BreachSet.SetParam gets 'combine' option
- BreachSamplesPlot improvement  
- Better support for req_monitors which are not stl_monitors
- Support for requirement parameters sampling at BreachRequirement level
- Updated Automatic transmission example

## Bugfixes
- Sim() recompute without argument now checks if traces are present before 
  recomputing with default time (was causing weird GetRobustSat behavior)
- log_format deactivated when using diskcaching wo parallel until figuring better solution 

## Developpers
- Added jsonlab toolbox 


# Release 1.3.3

- Fixed violation signals for alw (A => ev_[] B) types of requirements
- BrDemo are now all working with BreachRequirement class
- BreachProblem supports pre-conditions of BreachRequirement, treated
  as hard constraint for optimization
- Simulation is correctly skipped when pre-conditions on inputs only are not
  satisfied 
- overshoot requirement example added as alw (A=> ev_[] B) example for AFC
  
# Release 1.3.2

- plot_diagnosis now considers alw (A => ev_[] B) types of formulas 

# Release 1.3.1

- initial support of preconditions in BreachRequirement 
- stl_monitor plot_diagnosis method calls SplotSat
- BreachSamplesPlot plots sum when all true, all false or only one req

## Bugfixes

- bugfix for SaveInCache Falsification_runs for second and subsequent runs
- bugfix for index when computing signal signature with aliases

# Release 1.3.0

- BreachSamplesPlot plots bars and number of satisfied and violated
  requirement
- Alias (signal mapping) support fixed and improved
- Performance improvement for ExportTraces/SaveResults
- stl_monitor now should provide quant_sat and violation signals
  implemented for alw_monitor and alw_A_implies_B


# Release 1.3.beta4 
- BreachSamplesPlot provides interactive plots of BreachRequirement
  results
- BreachProblem saves runs in cache folders when DiskCaching is
  enabled
- runs of falsification problems can be loaded using
  FalsificationProblem.load_runs static method  
- bugfix for alw A=>B template and plot_diagnosis now separates A and
  A=>B on two axes


# Release 1.3.beta3
- DiskCaching for BreachProblem now save runs in the cache folders (Issue #29)
- skip computation of previously computed objective values in a BreachProblem in serial mode 
- fixed bug in BreachProblem.ResetObjective 
- alias appear in list in context menu BreachSignalsPlot
- Printing for requirement streamlined
- Printing signals now displays attributes and aliases 
- initial implementation of preconditions in BreachRequirement
- refactoring of BreachRequirement evalAllTraces, eliminate duplicates of input signals 
- support for multiple signal mappings/aliases 

# Release 1.3.beta2
- Renamed ogs into postprocess_signal_gen and formulas into
  req_monitors in BreachRequirement class
- STL_Eval more robust to NaN and Inf values 
- SplotSat now default to stairs plots until specified otherwise as
  global option


## Bugfixes
- STL_Formula copy bug resulting in stl_monitor not working 
- Display bug that was showing Nobody as name for BreachSimulinkSystem sometimes
- stl_monitor.eval not returning the correct value in X
- parameter lost from STL formula to stl_monitor


# Release 1.3.beta1
- New class BreachRequirement with hybrid and generic requirement
  evaluation support
- New plotting class BreachSignalsPlot
- Improvements in solver support for BreachProblem


# Release 1.2.18

- STL_Eval more robust to NaN and Inf values 
- SplotSat now default to stairs plots until specified otherwise as
  global option, also run InitBreach

# Release 1.2.17

## Bugfixes
- Fixed fix for STL monitoring interpolation bug (Issue #62) 

# Release 1.2.16
## Users
- printStatus now prints error message obtained during failed Simulink
  simulation
- New SimInModelsDataFolder option, default to false (i.e., by
  default, Breach won't change folder for simulalink simulations
- Changed default max_obj_eval to 100 and max_time to inf in BreachProblem

## Bugfixes 
- GetBrSet_Logged, _False, _True methods now works as expected in
  serial and parallel (BreachTEMA Issue #33)
- DiskCaching matfiles writable option should improve stability
- Fixed interpolation bug in STL_EvalThom (Issue #62)
- Fixed bug in BreachTraceSystem with useless legacy Xf field
- ModelsData folders properly created at initialization if not
  existant (BreachTEMA Issue #30)

## Developpers
- log_traces field for BreachProblem now consistent with serial and
  parallel (BreachTEMA Issue #31)
- StoreTracesOnDisk true by default for DiskCaching (BreachTEMA Issue #21)

# Release 1.2.15

## Users 
- New function FilterTracesStatus to remove/extract Null traces
  resulting from Simulink error and NaN traces resulting from input
  constraints violations 
- Simulink generated files are now all in ModelsData folder
- Added a small vdp example 
- In parallel, simulation are run in a different folder for each
  worker
- SaveResults now supports DiskCaching 


## Developpers
- In BreachProblem, objective function is now same for parallel and
  serial 

# Release 1.2.14

## Users
- DiskCaching feature added, enabled by SetupDiskCaching method
- Added BreachAbout, gives some info about which folder Breach runs
  from and its version (will add more in future releases)

## Bugfixes
- display status during optimization no longer display twice in
  certain cases
- parameters in workspace detected even if not part of Simulink model

## Developpers
- display_status_header is no longer called by FevalInit
- when using parallel computation, DiskCaching is used to log traces 


# Release 1.2.13

## Users
- Improves parallel support: display simulations during computation
  and optimization and number of workers now tunable with B.SetupParallel(NumWorkers)
- Experimental struct parameter detection now disabled by default
  (would not work with nested struct), made available with FindStruct
  option in BreachSimulinkSystem constructor 
- New signal generators (sinusoid, exponential, spike)

## Bugfixes
- freq_update now only affects display, independant from number of
  cores during parallel simulation

# Release 1.2.12 

## Users	
- InitBreach now better removes path from older/other versions
- InitBreach takes a second Boolean argument, if true, forces re-initialization
- Added BreachVersion, returns current version
- SaveResults, ExportToStruct, ExportToExcel now available for BreachSet and all derivated classes

## Bugfixes
- InitBreach creates ModelsData folder when absent, e.g., with fresh
  clone 

# Release 1.2.11

## Developpers
- ResetTimeSpent method for BreachProblem

## Bugfixes
- Problems when exporting to Excel with imported data (inputs or not, issues 56-57)
- BreachGui wouldn't start without sets in the workspace (issue 55) 
- Error message caught when checking with no trace 

# Release 1.2.11beta3
- Performance improvement for Export to Excel

## Developpers
- Some cleaning in BreachGui

## Bugfixes
- plot in GUI would not update axes labels and limits
- Fixed errors when clicking buttons and stuff on empty GUI
- Ignore file_idx parameter when writing traces 
- fixed bug in SampleDomain when combining samples of dimensions more than 1
- fixed bug in PlotDomain for more than three dimensions
- fixed load parameter set not working in GUI

# Release 1.2.11beta2

## Bugfixes
- small bug in displaying number of computed traces
- fixed bug in until robustness computation causing crash


# Release 1.2.11beta

## Users
- SetInitFn method sets a callback function for initialization before simulation
- new helper function isSignal returns if a name represents a signal
  in a BreachSet
- GUI safeguards around selected parameters for sampling (Issues 32-33-34)
- GUI 'Domain' button changed to 'Variables' 
- Support struct parameters and model workspace, still experimental 
- Enable parallel from GUI
- SaveSignals/Load
- BreachSave function saves Breach objects in the workspace into a mat file 
- GetVariables method returns names of variables in a set, i.e., parameters with non-empty domain/range

## Bugfixes
- GUI: fixed Env. Param button when no input signal
- Calling InitBreach in InstallBreach at very beginning to avoid
  missing varargin2struc (#52)
- Fixed bug in BreachDomain with 'int' and empty domain leading to change to enum type
- Fixed bug in SetInputGen which would not update domains appropriately

## Developpers
- New class BreachOptionGui to create gui from options
- FalsificationProblem now ignores requirement variables 
- ParamSynthesisProblem now ignores non-requirement variables

# Release 1.2.10

## Users
- InstallBreach option linear_interp, now default to 0 
- PrintAll and PrintParams now displays enum values
- Simulink Breach menu for creating interface improvements 

## Bugfixes
- GUI: selecting worst samples, etc, does not change plot parameters any more (#26-1)
- RobustAnd can't return inconsistent time values (not time advancing)
- GUI: embarrassing st_info bug when changing sampling option fixed (#26-4)

## Developpers
- PlotParams not using DiscrimPropValues any more
- More robust implementation of BreachDomain constructor 

# Release 1.2.9

## Users
- STL_ReadFile now supports parameters overriding existing functions, mex-files, etc (e.g., x[t]>eps) 
- PlotSatParams now shows superposed samples with different satisfaction values
- GUI can be opened by BreachGUI command, without argument
- GUI Check formula button sorts results by satisfaction and robustness
- PlotSatParams has datacursormode on by default, shows parameter values and robustness 
- SaveSignals method saves signals in a format that from_files_signal_gen can read easily 
- Many changes in GUI

## Bugfixes
- GetSignalsValues with itraj would assume itraj=1 
- ParamSynth/ReqMining works when changing the solver of param synth problem (Issue #10)

## Developpers
- GetSatValues returns unique values and params as well as all values and non-unique params
- BreachPlot and BreachPlotSat, new classes for listening plots

# Release 1.2.8

## Users
- Minor changes in Excel exported from SaveResults  
- GenDoc function to generate html documentation from a set of script

## Bugfixes
 - SaveResults now uses v7.3 format for large files
- Fixed issue in Autotrans_tutorial  

## Developpers
- new function get_breach_path 
- CleanModelsData erase stuff in Ext/ModelsData folder 

# Release 1.2.7

## Bugfixes
- GetSatValues not returning the proper values (caused issue #14)
- BreachProblem with 4 arguments not properly setting  domains (issue #19)
- Fixed issue with signal builder parameters that should be ignored in sim_breach (issue #19) 
- Fixed bug in BreachProblem that was converting singular domains into unbounded ones (issue #18)

## Developpers
- FEvalInit now updates obj_best and x_best fields rather than setting it

# Release 1.2.6

## Users 
- GetSignalValues can collect values from one or a subset of traces instead of all available

## Bugfixes
- Fixed ExportTracesToStruct  output (time series for data)
- Fixed Excel template path 
- Fixed PlotParam

# Release 1.2.5

## Bugfixes
- Fixed bug in STL_ExtractSignal (quotes around varname)
- Fixed PlotRobusMap changing order of props_names (Issue 14)
 
# Release 1.2.4

## Users 
- Changed CHANGELOG to CHANGELOG.md
- ExportTraces method for Simulink 
- SaveResults method creates a folder with  breach_system and summary files  and traces subfolder with trace files
- LoadResults reads a folder created with SaveResults

## Bugfixes
- Fixed bug in ComputeTraj (Issue 13)
- STL formulas now can use custom p-file, m-file, mex and builtin, though parameters with same name get shadowed   
- Issue with z-axis label in PlotSatParams
- GUI creating new set updates display
- Traj GUI: fixed bug preventing recomputation of trajectory when modifying a parameter 

## Developpers
- get_checksum method computes the checksum of Simulink mdl, and returned boolean true if it changed wrt creation
- mdl BreachSimulinkSystem property is a struct with model name, original path and date of creation
- varargin2struct function can be used to setup basic optional argument with syntax optionName, optionValue

# Release 1.2.3

## Users
- New method PlotSatParams to implement red-green dot maps available in GUI
- New method GetSetValues to get satifaction values of a previously monitored requirement

## Bugfixes
- fixed plot satisfaction map in GUI
- fixed PlotRobustMap method

# Release 1.2.2

## Users
- ReqMiningProblem constructor with two argument (synth_pb and falsif_pb, e.g.) now works no matter the order of argument, but will error if it cannot identify a synthesis problem and a falsification problem
- SetupLogFolder default now has date of creation in name in format mdl_name_ddmmyy-HHMM 
- Predicates accept builtin and custom  user defined functions 

## Bugfixes
- fixed problem when shift-selecting multiple parameters in main GUI  
- fixed falsification (and all types of problems) not enforcing domains properly 
- fixed renaming in GUI not renaming in base workspace
- fixed GUI allows multiple cell edition with shift-return

## Developpers
- BreachProblem now stores domains for params  
- BreachProblem has a CheckinDomain method 

# Release 1.2.1

## Users
- GUI: set input now updates main GUI when done
- InitNNDemo in simple ML model example
- GUI allows multiple cell edition with shift-return
- GUI: set input does not reset signal generators

## Bugfixes
- GUI copy select button works again

## Developpers

# Release 1.2.0

## Users
- New PlotDomain method
- domains can be set with the syntax SetDomain([4 5.1]), will infer double type
- error is returned when wrong type is given in SetDomain
- default behavior of legacy sampling functions is to reset and replace, this is changed by AppendWhenSample flag
- new warning when breach can't get checksum of model from Simulink. Only affects the Log to file features.
- GetBoundedDomains returns domains with bounds/ranges
- GUI work with all BreachSet in the workspace, regardles of system interfaced
- AddSpec now add property parameters automatically
- PlotParams, PlotSignals, etc, now re-use the current figure if there is one (gcf)

## Bugs fixed
- compilation issue with VS 2015 for online monitoring (added include minmax and algorithm) 
- from_file_signal_gen not working for files in current folder 
- AppendWhenSample now works with legacy sampling functions, i.e., concat new param vectors when true 
- Parallel with less initial parameters than freq_update now works as intended
- AddSpec now checks whether formula is already present (prevent overwrite param values)

## Developpers  
- SampleDomain now calls sample method of BreachDomain
- BreachSet now has only one field Domains for signal and parameter domains 
- SetParamRanges errors when param does not exist (before would create it)
- ResetParamSet uses Domains instead of ParamRanges field
- SetParamSpec has an additional argument to ignore if param is a system parameter, if not set, return an error
- got rid of ParamRanges field, useless since we use domains

# Release 1.1.0

## Users 
- New menu in Simulink editor, can generate a Breach interface
- Breach now supports/detects To Workspace blocks
- Breach now supports/detects signal builders
- GUI for specs now can change parameters of formulas
- New GUI for creating input generator, accessible from main GUI 
- Main GUI keyboard shortcuts
- rightarrow from constant parameters to varying parameters, equivalent to [add =>]
- leftarrow from varying parameters equivalent to button [<= rem]
- return from varying parameters focus on min value
- log to files via SetupLogFolder 
- Support for different types via Domains
- SampleDomain, more comprehensive sampling of domains

## BugFix
- Fixed problem when resetting input generator with models with no inputs
- Fixed problem when logging multi dimensional signals
- Fixed problem with empty tspan when interfacing model
- Fixed problem Simulink models with no inputs


# Release 1.0.0

## Users

- BrDemo subfolder Optim contains InitAFC_Falsif which creates falsification problems
- BrDemo has a subfolder InputGen containing demo script for each type of generator (under construction)
- New method FilterSpec return sets that satisfy and don't satisfy a spec  
- SetParamRanges accepts a single interval for multiple parameters, 
    Example:
B.SetParamRanges({'p1', 'p2','p3'}, [ -1 , 1] ) is equivalent to
B.SetParamRanges({'p1', 'p2','p3'}, [ -1 , 1 ; -1 1 ;-1 1] ) 

- BreachSimulinkSystem now has a field sim_args used to pass additional option arguments to the Simulink sim command
- BreachSystem has SetSpec method - same as AddSpec, but reset Specs before
- BreachOpenSystem (hence BreachSimulinkSystem) has a AddInputSpec and SetInputSpec, calls AddSpec or SetSpec its InputGenerator
- PrintSignals now says which signal is an input 
- GetInputSignalsIdx returns indices of input signals in the list of signals of a system
- Minor improvement on display status for BreachProblem
- Input specs now are taken into account in Falsification problems
- Breach now stores its model copy and GUI files in a unique folder in Ext/ModelsData
