t = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
v = [-2, -3, -4, 1, 2, 3, -1, -1, -1, 1, 2, 3, 4];

s.times = t;
s.values = v;

interval.begin = 1;
interval.end = 4;

bound.begin = 1;
bound.end = 3;

implicant = BreachImplicant;
implicant = implicant.addInterval(interval.begin, interval.end);

target = 2;

[out_implicant, error] = BreachDiagnostics.diag_ev_t(s, implicant, bound, target);

%signal.times = t;
%signal.values = v;

%interval.begin = 2.5;
%interval.end = 4;

%target_value = 2;

%t2 = [1.5, 2.5];

%implicant = BreachImplicant;
%implicant = implicant.addInterval(3,2);
%implicant = implicant.addInterval(3,5);
%implicant = implicant.addInterval(-inf,-3);
%implicant = implicant.setWorstTime(3.4);
%implicant.getInterval(1)
