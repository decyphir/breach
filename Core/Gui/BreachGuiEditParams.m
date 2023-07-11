classdef BreachGuiEditParams < BreachGuiClass
 % Defines a Gui class that takes a BreachSet as argument and implements:
 %  - a listbox_params 
 %  - a slider_samples
 %  - a panel_params 
 %  - a slider_edit_params connected to the above
    methods
        function this =BreachGuiEditParams(B)

            if nargin==0
                % default BreachSet for testing purposes
                B = BreachSet({'p1','p2', 'Some_other_param'});
                B.SetParamRanges('p1', [0  12]);
                B.SetParamRanges('p2', [-1 5]);
                B.QuasiRandomSample(10);
            end
            %
            
            this.data_gui.BrSet =B;
            this.data_gui.all_params = B.GetParamList();
            this.data_gui.current_sample = 1;             
           
            % listbox
            id_listbox = 'listbox_params';
            string = this.get_string_params_and_values(this.data_gui.all_params, this.data_gui.current_sample);
            c = @(o,e)(this.callback_listbox(id_listbox));
            this.create_listbox(id_listbox,string, c);
            
            
            % slider for sample idx
            max_slider= B.GetNumSamples();
            %if max_slider>1
            min_slider = 1;
            
            id_slider_samples= 'slider_samples' ;
            c_slider_samples = @(o,e)(this.callback_slider_samples());
            e_slider =            this.create_slider(id_slider_samples, 'what', c_slider_samples);
            e_slider.h = 0.5;
                        
            set(e_slider.hdle,'Min', min_slider);
            set(e_slider.hdle,'Max', max_slider);
            set(e_slider.hdle,'Value', 1);
            set(e_slider.hdle,'SliderStep',[1 10]/(min(max_slider, min_slider+1) - min_slider));
            
            panel_layout = { {id_slider_samples}; {id_listbox}};
            %else
           %     panel_layout = {{id_listbox}};
           % end
            % embed into a panel (nicer)
            this.create_panel('panel_params', 'Parameters',panel_layout);
            
            this.update_panel_params();
            % slider 
            id_slider = 'slider_edit_param';
            name = this.data_gui.all_params{1};
            domain = B.GetDomain(name);
            value=  B.GetParam(name, 1);
            c_slider = @this.callback_slider;   
            
            this.create_domain_slider(id_slider, name,  domain,value,c_slider);
            
            % layout
            layout = {{'panel_params'} ;...
                {id_slider} };

            % update
            this.update_domain_slider('slider_edit_param', false);

            this.set_layout(layout);
            this.enable_resizable();
        end

        function callback_listbox(this,id)
            e = this.uimap(id); % listbox elem                
            idx_param = get(e.hdle,'Value');
            param = this.data_gui.all_params{idx_param};

            domain  = this.data_gui.BrSet.GetDomain(param);

            this.data_gui.slider_edit_param.name = param;      
            this.data_gui.slider_edit_param.domain = domain;
            this.data_gui.slider_edit_param.value = this.data_gui.BrSet.GetParam(param,this.data_gui.current_sample);
            this.update_domain_slider('slider_edit_param', false);
        end
        
        
        function callback_slider_samples(this)
            e = this.uimap('slider_samples');
            val = round(get(e.hdle, 'Value'));
            this.data_gui.current_sample=val;            
            set(e.hdle,'Value', val);
            % 
            e_listbox = this.uimap('listbox_params'); % listbox elem                
            idx_param = get(e_listbox.hdle,'Value');
            param = this.data_gui.all_params{idx_param};
            
            % update slider edit param
            this.data_gui.('slider_edit_param').value = this.data_gui.BrSet.GetParam(param,this.data_gui.current_sample);
            this.update_domain_slider('slider_edit_param', false);

            this.update_panel_params();            
            
            
        end


        
        function callback_slider(this)
            e = this.uimap('slider_edit_param_slider');

            val = get(e.hdle, 'Value');
            param = this.data_gui.slider_edit_param.name;
            all_params_values = this.data_gui.BrSet.GetParam(param);
            all_params_values(this.data_gui.current_sample) = val; 
            this.data_gui.BrSet.SetParam(param, all_params_values);
            this.update_panel_params();

        end


        function update_panel_params(this)
            e = this.uimap('panel_params');            
            title_panel = ['Parameters for sample ' num2str(this.data_gui.current_sample) '/' num2str(this.data_gui.BrSet.GetNumSamples())];
            set(e.hdle, 'Title', title_panel);

            string_listbox = this.get_string_params_and_values(this.data_gui.all_params, this.data_gui.current_sample);
            e_listbox = this.uimap('listbox_params');
            set(e_listbox.hdle,'String', string_listbox);

        end

        function update_samples_slider(this)
            val = this.data_gui.current_sample;
            num_samples =this.data_gui.BrSet.GetNumSamples();
        end        
         
        
        function str = get_string_params_and_values(this,params, idx_sample)
            
            str = cell(1, numel(params));
            for idx_param = 1:numel(params)
                param = params{idx_param};
                val = this.data_gui.BrSet.GetParam(param, idx_sample);
                str{idx_param}=  [param  ': ' num2str(val)];
            end
        end

    end
end
