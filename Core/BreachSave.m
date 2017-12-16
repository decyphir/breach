function tosave = BreachSave(filepath)
% BreachSave(filepath) saves workspace Breach objects

var = evalin('base', 'who');
tosave = {};
for iv = 1:numel(var)
    vari = var{iv};
    if isa_breach_obj(vari)
        tosave =union(tosave, ['''' vari '''']);
    end
end

cmd = ['save(''' filepath ''',' strjoin(tosave,',') ');'];
evalin('base', cmd);
    
end

function b = isa_breach_obj(v)

c = evalin('base', ['class(' v ')' ]);
b = isequal(c, 'STL_Formula')||...
    ~isempty(regexp(c, '^Breach\w+', 'once'))||...
    evalin('base',['isa(' v ',''BreachProblem'')'])||...
    evalin('base',['isa(' v ',''BreachStatus'')']);

end