% STUDENT DETAILS
%{
                ENM221-0055/2017
                Gladys Wanjeri Gachoka                
                B.Sc. Mechatronic Engineering 
                EMT 2436 : Vibrations
                
                11 June, 2021
%}

% LAB WORK 1
%{                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Create a program to solve second order Constant Coefficient Homogenous ODEs 
(without using dsolve). Ask user for the 3 contants (a, b, c) and the 2
initial conditions.
%}

clc

syms y(t) C1 C2 C3 C4 C5 C6

% Asking for the 3 constant coefficients & 2 input conditions & storing
% them in variables
disp("Constant Coefficient Homogenous ODE");
coeff_a = input('Enter coeeficient a : ');
coeff_b = input('Enter coeeficient b : ');
coeff_c = input('Enter coeeficient c : ');
fprintf("\n");
disp("Assuming all initial conditions are at t = 0");
cond_1 = input('Enter condition 1 : y(0) = ');
cond_2 = input('Enter condition 2 : y''(0) = ');
fprintf("\n");

% Creating the ODE
ode = coeff_a * diff(y,t,2) + coeff_b * diff(y,t) + coeff_c * y == 0;

% Determining nature of the roots hence nature of general sln
% If b^2 - 4ac > 0 ; Real roots
% If b^2 - 4ac = 0 ; Real repeated roots
% If b^2 - 4ac < 0 ; Complex roots
root_type = coeff_b^2 - 4 * coeff_a * coeff_c;

if ~isempty(cond_1) && ~isempty(cond_2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% REAL ROOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if root_type > 0
        disp("The equation has real roots");
        roots = quadratic(coeff_a,coeff_b,coeff_c);
        disp("Roots : ");
        disp(roots);

        % Real general solution
        genEqn1 = C1 * exp(roots(1)*t) + C2 * exp(roots(2)*t);
        disp("General equation : ");
        disp(genEqn1);

        % Applying initial conditions    
        c1 = C1 * exp(roots(1)*t) + C2 * exp(roots(2)*t) == cond_1; % applying condition 1
        s1 = subs(c1, t, 0); % expr when t = 0
        disp("Eqn 1 : ");
        disp(s1);
        cc1 = isolate(s1, C1); % making C1 the subject

        c2 = diff(genEqn1) == cond_2; % applying condition 2
        s2 = subs(c2, t, 0); % expr when t = 0
        disp("Eqn 2 : ");
        disp(s2);
        cc2 = subs(s2, lhs(cc1), rhs(cc1)); % substituting C1 with its value in 
                                           % cc1 to get eqn with C2 as only unknown

        valC2 = solve(cc2, C2, 'Real', true); % solving cc2 to get value of C2
        disp('C2 = ');
        disp(valC2);
        s3 = subs(s1, C2, valC2); % substituting C2 with its real value valC2 
                                 % in s1
        valC1 = solve(s3, C1, 'Real', true); % solving s3 to get value of C1
        disp("C1 = ");
        disp(valC1);

        % Final answer
        sln1 = valC1 * exp(roots(1)*t) + valC2 * exp(roots(2)*t);
        disp("Answer is : ");
        disp(sln1);

        % NB : vpa(x, d) when numbers get too large/strange

    %%%%%%%%%%%%%%%%%%%%%%%%%% REAL REPEATED ROOTS %%%%%%%%%%%%%%%%%%%%%%%%
    elseif root_type == 0 
        disp("The equation has real repeated roots");
        roots_2 = quadratic(coeff_a,coeff_b,coeff_c);
        disp("Roots : ");
        disp(roots_2);

        % Real repeated general solution
        genEqn2 = C3 * exp(roots_2(1)*t) + C4 * t * exp(roots_2(2)*t);
        disp("General equation : ");
        disp(genEqn2);

        % Applying initial conditions
        c1_2 = C3 * exp(roots_2(1)*t) + C4 * t * exp(roots_2(2)*t) == cond_1; % applying condition 1
        s1_2 = subs(c1_2, t, 0); % expr when t = 0
        disp("Eqn 1 : ");
        disp(s1_2);
        cc1_2 = isolate(s1_2, C3); % making C3 the subject

        c2_2 = diff(genEqn2) == cond_2; % applying condition 2
        s2_2 = subs(c2_2, t, 0); % expr when t = 0
        disp("Eqn 2 : ");
        disp(s2_2);
        cc2_2 = subs(s2_2, lhs(cc1_2), rhs(cc1_2)); % substituting C3 with its value in 
                                                   % cc1_2 to get eqn with C4 as only unknown

        valC4 = solve(cc2_2, C4, 'Real', true); % solving cc2_2 to get value of C4
        disp('C4 = ');
        disp(valC4);
        s3_2 = subs(s1_2, C4, valC4); % substituting C4 with its real value valC4
                                     % in s1_2
        valC3 = solve(s3_2, C3, 'Real', true); % solving s3_2 to get value of C3
        disp("C3 = ");
        disp(valC3);

        % Final answer
        sln2 = valC3 * exp(roots_2(1)*t) + valC4 * t * exp(roots_2(2)*t);
        disp("Answer is : ");
        disp(sln2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLEX ROOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif root_type < 0
        disp("The equation has complex roots");
        roots_3 = quadratic(coeff_a,coeff_b,coeff_c);
        disp("Roots : ");
        disp(roots_3);

        n = real(roots_3(1));
        m = imag(roots_3(1));

        % Complex general solution
        genEqn3 = exp(n*t) * (C5 * cos(m*t) + C6 * sin(m*t));
        disp("General equation : ");
        disp(genEqn3);

        % Applying initial conditions
        c1_3 = exp(n*t) * (C5 * cos(m*t) + C6 * sin(m*t)) == cond_1; % applying condition 1
        s1_3 = subs(c1_3, t, 0); % expr when t = 0
        disp("Eqn 1 : ");
        disp(s1_3);
        cc1_3 = isolate(s1_3, C5); % making C5 the subject

        c2_3 = diff(genEqn3) == cond_2; % applying condition 2
        s2_3 = subs(c2_3, t, 0); % expr when t = 0
        disp("Eqn 2 : ");
        disp(s2_3);
        cc2_3 = subs(s2_3, lhs(cc1_3), rhs(cc1_3)); % substituting C5 with its value in 
                                                   % cc1_3 to get eqn with C6 as only unknown

        valC6 = solve(cc2_3, C6, 'Real', true); % solving cc2_3 to get value of C6
        disp('C6 = ');
        disp(valC6);
        s3_3 = subs(s1_3, C6, valC6); % substituting C6 with its real value valC6
                                     % in s1_3
        valC5 = solve(s3_3, C5, 'Real', true); % solving s3_3 to get value of C5
        disp("C5 = ");
        disp(valC5);

        % Final answer
        sln3 = exp(n*t) * (valC5 * cos(m*t) + valC6 * sin(m*t));
        %sln3b = exp(n*t) * valC5 * cos(m*t) + exp(n*t) * valC6 * sin(m*t); % alternative
        disp("Answer is : ");
        disp(sln3);  

    else
        disp("One or more coefficients missing!");

    end
else
    disp("One or more initial conditions missing!");
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHORTCUT SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(root_type) && ~isempty(cond_1) && ~isempty(cond_2)
    Dy = diff(y);
    cond1 = y(0) == cond_1;
    cond2 = Dy(0) == cond_2;
    conds = [cond1 cond2];

    sln = dsolve(ode, conds);

    disp("Sln using dsolve = ");
    disp(sln);
end

























