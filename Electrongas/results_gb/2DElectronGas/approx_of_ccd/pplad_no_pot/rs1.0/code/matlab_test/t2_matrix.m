function t2mat = t2_matrix(t2mat_old)

%
% Set up a new coupled-cluster T2 matrix given an old matrix
% t2mat_old.
%

% Interaction matrix elements
v0125 = 0.313328534432575;
v0134 = -0.313328534432575;
v0101 = 1.25331413773030;
v2525 = 0.861653469911781;
v2534 = -0.234996398659501;
v3434 = 0.861653469911781;
v1515 = 0.626657068865150;
v0505 = 0.939985603297726;
v0202 = 0.626657068865150;
v1212 = 0.939985603297726;
v1414 = 0.939985603297726;
v0404 = 0.626657068865150;
v0303 = 0.939985603297726;
v1313 = 0.626657068865150;

% Fock matrix elements
f00 = 1.0 + 1.25331413773030;
f11 = 1.0 + 1.25331413773030;
f22 = 2.0 + 0.626657068865150 + 0.939985603297726;
f33 = 2.0 + 0.939985603297726 + 0.626657068865150;
f44 = 2.0 + 0.626657068865150 + 0.939985603297726;
f55 = 2.0 + 0.939985603297726 + 0.626657068865150;

% Energy denominator
denom2501 = f00 + f11 - f22 - f55;
denom3401 = f00 + f11 - f33 - f44;

t2501_old = t2mat_old(1);
t3401_old = t2mat_old(2);

%----------------------------------------------------------------
%
% t2501
%
%----------------------------------------------------------------

% Diagram 1
diagram1(1) = v0125;

% Diagram 2
diagram2(1) = v2525*t2501_old + v2534*t3401_old;

% Diagram 3
diagram3(1) = (v0101 + v0125*t2501_old + v0134*t3401_old)*t2501_old;

% Diagram 4
diagram4(1) = (-v1515 + 0.5*v0125*t2501_old)*t2501_old ...
    + (-v0505 + 0.5*v0125*t2501_old)*t2501_old ...
    + (-v0202 + 0.5*v0125*t2501_old)*t2501_old ...
    + (-v1212 + 0.5*v0125*t2501_old)*t2501_old;

% Diagram 5
diagram5(1) = -v0125*t2501_old*t2501_old ...
        - v0125*t2501_old*t2501_old

% Diagram 6
diagram6(1) = -2.0*v0125*t2501_old*t2501_old ...
        - 2.0*v0134*t3401_old*t2501_old;

%----------------------------------------------------------------
%
% t3401
%
%----------------------------------------------------------------

% Diagram 1 
diagram1(2) = v0134;

% Diagram 2
diagram2(2) = v2534*t2501_old + v3434*t3401_old;

% Diagram 3
diagram3(2) = (v0101 + v0125*t2501_old + v0134*t3401_old)*t3401_old;

% Diagram 4
diagram4(2) = (-v1414 + 0.5*v0134*t3401_old)*t3401_old ...
        + (-v0404 + 0.5*v0134*t3401_old)*t3401_old ...
        + (-v0303 + 0.5*v0134*t3401_old)*t3401_old ...
        + (-v1313 + 0.5*v0134*t3401_old)*t3401_old;

% Diagram 5
diagram5(2) = -v0134*t3401_old*t3401_old ...
        - v0134*t3401_old*t3401_old

% Diagram 6    
diagram6(2) = -2.0*v0125*t2501_old*t3401_old ...
        - 2.0*v0134*t3401_old*t3401_old;
    
   
%----------------------------------------------------------------

temp = diagram1(1) + diagram2(1) + diagram3(1) + diagram4(1) ...
        + diagram5(1) + diagram6(1);
t2mat(1) = temp/denom2501;

temp = diagram1(2) + diagram2(2) + diagram3(2) + diagram4(2) ...
        + diagram5(2) + diagram6(2);
t2mat(2) = temp/denom3401;
