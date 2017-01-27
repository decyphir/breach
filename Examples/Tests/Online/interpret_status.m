function st = interpret_status(report)

status = report.status;
r = report.online.rob.values(1);
rob_thom = report.thom.values(1);
st = '';
if status(1)==0
    st=sprintf('!! Robustness online %g and Thom %g are different (diff:%g).\n', r(1), rob_thom(1),r(1)-rob_thom(1));
end
if status(2)==0
    st=sprintf('!! Lower robustness greater than upper.\n');
end
if status(3)==0
    st=sprintf('!! Lower robustness greater than rob.\n');
end
if status(4)==0
    st=sprintf('!! Robustness greater than upper rob\n');
end
if status(5)==0
    st=sprintf('!! Upper rob not decreasing\n');
end
if status(6)==0
    st=sprintf('!! lower rob not increasing\n');
end
if all(status)
    st=sprintf('ALL GOOD.\n');
end

end
