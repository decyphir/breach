function signals = STL_ExtractSignals(phi)
%STL_EXTRACTSignals extracts names of signals  involved in phi
%
% Synopsis: signals = STL_ExtractSignals(phi)
%
% Input:
%  - phi : the formula from which the signals are extraced
%
% Output:
%  - signals: list of signals involved in phi (note: currently only detec-
%    ted from patterns of the form signal_id[t])
%

[~, signals] = STL_ExtractPredicates(phi);