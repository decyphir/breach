classdef postprocess_fcn_wrapper < output_gen
    properties
        fun
    end
    
    methods
        function this = postprocess_fcn_wrapper(fn_handle,  sigs_in, sigs_out, params, p0)
            this.fun = fn_handle;
            this.signals_in = sigs_in;
            this.signals = sigs_out;
            if exist('params','var')&&~isempty(params)
                this.params = params;
                this.p0 = p0;
            end
        end
        
        function [time, X]= computeSignals(this, time,Xin, p)
            
            eval_st = '[';
            for ia = 1:numel(this.signals)-1
                eval_st  = [ eval_st 'X(' num2str(ia) ',:),'];
            end
            eval_st = [eval_st 'X(' num2str( numel(this.signals)) ',:)] = this.fun(time,' ];
            %    eval_st = [eval_st 'X(' num2str( numel(this.signals)) ',:)] = this.fun(' ];
            
            for ia = 1:numel(this.signals_in)-1
                eval_st  = [ eval_st 'Xin(' num2str(ia) ',:),'];
            end
            eval_st = [eval_st 'Xin(' num2str( numel(this.signals_in)) ',:)'];
            
            for ip = 1:numel(p)
                eval_st = [eval_st ',p(' num2str(ip) ')'];
            end
            eval_st = [eval_st ');'];
            eval(eval_st);
        end
        
        
    end
end