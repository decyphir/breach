function trace = breach_result_rename(trace, varargin)

        arg_err_msg = 'Argument should be a list of pairs of signal names, or two cells of signal names with same size.';
            switch nargin
                case 3
                    if iscell(varargin{2})
                        if ~iscell(varargin{2})||(numel(varargin{1}) ~= numel(varargin{2}))
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                        for is = 1:numel(varargin{2})
                            rename(varargin{1}{is},  varargin{2}{is});
                        end
                    else
                        if ischar(varargin{1})&&ischar(varargin{2})
                            rename(varargin{1}, varargin{2});
                        else
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
                otherwise
                    for is = 1:numel(varargin)/2
                        try
                            rename(varargin{2*is-1}, varargin{2*is});
                        catch
                            error('SetSignalMap:wrong_arg', arg_err_msg);
                        end
                    end
            end
   
            
    function rename(sig_old, sig_new)
        idx =strcmp(sig_old, trace.outputs.names);
        if any(idx)
            trace.outputs.names{idx} = sig_new;
        end
        idx =strcmp(sig_old, trace.inputs.names);
        if any(idx)
            trace.inputs.names{idx} = sig_new;
        end
        
    end
end