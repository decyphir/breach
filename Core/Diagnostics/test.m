t1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
v1 = [-2, -3, -4, 1, 2, 3, -1, -1, -1, 1, 2, 3, 4];

t2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
v2 = [2, 2, 0, 1, -2, -2, 4, 2, 2, 2, -1, -1, -1];

s1.times = t1;
s1.values = v1;

s2.times = t2;
s2.values = v2;

interval.begin = 0;
interval.end = 2;

bound.begin = 1;
bound.end = 2;


implicant = BreachImplicant;
implicant = implicant.addInterval(interval.begin, interval.end);

target = 2;

[implicant1, implicant2] = BreachDiagnostics.diag_or_t(s1, s2, implicant, target);

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
