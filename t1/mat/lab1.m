close all
clear all

pkg load symbolic;

% minimum student number: 86361

R1 = sym ('1.0407324334365136');
R2 = sym ('2.0857867519728686');
R3 = sym ('3.071199615896663');
R4 = sym ('4.041770723234123');
R5 = sym ('3.140660073703873');
R6 = sym ('2.0936811661064763');
R7 = sym ('1.0443244500565343');
Va = sym ('5.194209863050843');
Id = sym ('1.0412616050274464');
Kb = sym ('7.021062278588699');
Kc = sym ('8.137326206873837');

% Mesh Method
syms IA IB IC ID

Eq_A = R1*IA + R3*(IA+IB) + R4*(IA+IC) == Va;
Eq_C = R6*IC + R7*IC + R4*(IC+IA)== Kc*IC;
Eq_B = IB == Kb*R3*(IA+IB);
Eq_D = ID == Id;

sol = solve(Eq_A,Eq_B,Eq_C,Eq_D);

% Sym to double
IA = double(sol.IA); 
IB = double(sol.IB);
IC = double(sol.IC);
ID = double(sol.ID);

Ib_m = IB;
Id = double(Id);

%Components current calculation
I_r1m = -IA;
I_r2m = IB;
I_r3m = IA + IB;
I_r4m = IA + IC;
I_r5m = ID - IB;
I_r6m = -IC;
I_r7m = -IC;

%Nodes voltage calculation (Ohm's Law)
V1_m = double(Va);
V2_m = V1_m - double(R1) * IA;
V3_m = V2_m + double(R2) * IB;
V4_m = double(R4) * (IA + IC);
V5_m = V4_m + double(R5) * (ID - IB);
V6_m = V4_m - double(Kc) * IC;
V7_m = - double(R6) * IC;
V8_m = V7_m;

% -----------------------------------------------------------------------------
%Nodal Method 
syms V0_n V1_n V2_n V3_n V4_n V5_n V6_n V7_n V8_n

Eq_1 = V8_n == V7_n;
Eq_2 = V0_n == 0;
Eq_3 = (V1_n - V2_n)/R1 +(V0_n - V8_n)/R6 + (V0_n - V4_n)/R4 == 0;
Eq_4 = (V2_n-V4_n)/R3 + (V2_n-V3_n)/R2 + (V2_n-V1_n)/R1 == 0;
Eq_5 = - Kb*(V2_n-V4_n) + (V3_n-V2_n)/R2 == 0;
Eq_6 = (V5_n-V4_n)/R5 + Kb*(V2_n-V4_n) - Id == 0;
Eq_7 = (V7_n-V6_n)/R7 + (V7_n - V0_n)/R6 == 0;
Eq_8 = V1_n - V0_n == Va;
Eq_9 = V4_n - V6_n == Kc * (V0_n - V7_n)/R6;

sol_n = solve(Eq_1,Eq_2,Eq_3,Eq_4,Eq_5,Eq_6,Eq_7,Eq_8,Eq_9);

% Sym to double
V0_n = double(sol_n.V0_n);
V1_n = double(sol_n.V1_n);
V2_n = double(sol_n.V2_n);
V3_n = double(sol_n.V3_n);
V4_n = double(sol_n.V4_n);
V5_n = double(sol_n.V5_n);
V6_n = double(sol_n.V6_n);
V7_n = double(sol_n.V7_n);
V8_n = V7_n;

I_r1n = (V2_n - V1_n)/double(R1);
I_r2n = (V3_n - V2_n)/double(R2);
I_r3n = (V2_n - V4_n)/double(R3);
I_r4n = (V4_n - V0_n)/double(R4);
I_r5n = (V5_n - V4_n)/double(R5);
I_r6n = (V8_n - V0_n)/double(R6);
I_r7n = (V6_n - V7_n)/double(R7);
Ib_n = I_r2n;

% -----------------------------------------------------------------------------
%Printing the values in a table

diary "final_table.tex"
diary on
printf("@Gb & %f & %f z\n", Ib_m, Ib_n);
printf("@id & %f & %f z\n", Id, Id);
printf("@r1 & %f & %f z\n", I_r1m, I_r1n);
printf("@r2 & %f & %f z\n", I_r2m, I_r2n);
printf("@r3 & %f & %f z\n", I_r3m, I_r3n);
printf("@r4 & %f & %f z\n", I_r4m, I_r4n);
printf("@r5 & %f & %f z\n", I_r5m, I_r5n);
printf("@r6 & %f & %f z\n", I_r6m, I_r6n);
printf("@r7 & %f & %f z\n", I_r7m, I_r7n);
printf("V1 & %f & %f z\n", V1_m, V1_n);
printf("V2 & %f & %f z\n", V2_m, V2_n);
printf("V3 & %f & %f z\n", V3_m, V3_n);
printf("V4 & %f & %f z\n", V4_m, V4_n);
printf("V5 & %f & %f z\n", V5_m, V5_n);
printf("V6 & %f & %f z\n", V6_m, V6_n);
printf("V7 & %f & %f z\n", V7_m, V7_n);
printf("V8 & %f & %f z\n", V8_m, V8_n);

diary off


