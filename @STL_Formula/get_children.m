function subphis = get_children(this)

subphis = [];
if ~isempty(this.phi) 
    subphis = {phi};
elseif ~isempty(this.phi1)
    subphis = {this.phi1, this.phi2};
end