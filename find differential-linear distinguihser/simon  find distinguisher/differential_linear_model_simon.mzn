include "Simon_diff_operations.mzn";
include "Simon_continuous_operations.mzn";
include "Simon_linear_operations.mzn";
%include "Simon_common_operations.mzn";



int: word_size = 16; 
% 差分部分的参数
int: nr_diff = 3;
array [0..nr_diff, 0..word_size-1] of var 0..1: left_diff;  % Left word of the state, input of the round
array [0..nr_diff, 0..word_size-1] of var 0..1: right_diff;  % Right word of the state, input of the round

array [0..nr_diff-1, 0..word_size-1] of var 0..1: varibits;  % variable used to store the sum of probability in the differential part 
array [0..nr_diff-1, 0..word_size-1] of var 0..1: doublebits;  % variable used to store the sum of probability for the differential part
array [0..nr_diff-1, 0..word_size-1] of var 0..1: Belt; % variable of output of SIMON_like_round_f
array [0..nr_diff-1, 0..word_size-1] of var 0..1: Z;
array [0..nr_diff-1, 0..word_size-1] of var 0..1: dummy_and;  % dummy variables for the RotationAnd operator in the differential part
array [0..nr_diff-1, 0..word_size-1] of var 0..1: CC;  % dummy variables for the RotationAnd operator in the differential part
var 0..4*word_size: diff_probability; % variable used to store the sum of probability

% 差分-线性部分的参数
int:nr_diffLin = 3; % Number of rounds for the differential-linear part
array [0..nr_diffLin, 0..word_size-1] of var lower..upper: diffLin_left;  % Left word of the state, at the beginning of the round for the differential-linear part
array [0..nr_diffLin, 0..word_size-1] of var lower..upper: diffLin_right;  % Right word of the state, at the beginning of the round for the differential-linear part
var lower..upper: diffLin_corr; % variable used to store the correlation for the middle part
var lower..upper: diffLinComplement_corr; % variable to store the approximation of f(r) = -log2(|r|) 

% 线性部分的参数
int: nr_lin = 3; % Number of rounds for the linear part (bottom part)
array [0..nr_lin, 0..1, 0..word_size-1]  of var 0..1: lin_state; 
array [0..nr_lin, 0..word_size-1]  of var 0..1: lin_and1; % dummy variable for the fork operation 
array [0..nr_lin, 0..word_size-1]  of var 0..1: lin_and2; % dummy variable for the fork operation  
array [0..nr_lin, 0..word_size-1] of var 0..1: p_lin; % variables to store the probability for the linear part
var 0..4*word_size: lin_corr; 


% 差分部分建模
predicate diff_round (
  array[0..word_size-1] of var 0..1: x0, 
  array[0..word_size-1] of var 0..1: x1, 
  array[0..word_size-1] of var 0..1: v, 
  array[0..word_size-1] of var 0..1: doub, 
  array[0..word_size-1] of var 0..1: b,
  array[0..word_size-1] of var 0..1: z,
  array[0..word_size-1] of var 0..1: d,
  array[0..word_size-1] of var 0..1: y0, 
  array[0..word_size-1] of var 0..1: y1,   
  int: r
) = (
  let {
    array [0..word_size-1] of var 0..1: x2 = LRot(x0, 8),
    array [0..word_size-1] of var 0..1: x3 = LRot(x0, 1),
    array [0..word_size-1] of var 0..1: x4 = LRot(x0, 2),
    array [0..word_size-1] of var 0..1: x5 = LRot(x0, 15),
    array [0..word_size-1] of var 0..1: x6 = LRot(x0, 7),
  } in
  SIMON_like_round_f(x0, v, doub, b, z) /\
  Xor_word(b, x1, y0, d, r) 
  
  
);

constraint forall (i in 1..nr_diff) (
 diff_round(
   row(left_diff,i-1), 
   row(right_diff,i-1), 
   row(varibits, i-1),
   row(doublebits, i-1),
   row(Belt, i-1),
   row(Z, i-1),
   row(dummy_and, i-1),
   row(left_diff,i), %注意左右交换
   row(right_diff,i), 
   i
 ) /\
 row(right_diff,i)=row(left_diff,i-1)
) /\ (sum((row(left_diff,0))) + sum(row(right_diff,0)) > 0)
/\ (sum((row(left_diff,0))) < word_size);

% 差分-线性部分建模
function array [int] of var lower..upper: continuous_round(array [int] of var lower..upper: xx, array [int] of var lower..upper: yy) = 
let { 
    int: n = length(xx);
    array [int] of var lower..upper: rr1 = continuous_LRot(xx, 1);
    array [int] of var lower..upper: rr8 = continuous_LRot(xx, 8);
    array [int] of var lower..upper: rr2 = continuous_LRot(xx, 2);
} in
array1d(
  0..2*word_size-1, 
  continuous_xor(continuous_xor(continuous_and(rr1, rr8), yy),rr2)++ xx   
);


constraint forall (i in 0..nr_diffLin-1) (
 continuous_round(row(diffLin_left,i),row(diffLin_right,i)) = 
   array1d(0..2*word_size-1, row(diffLin_left,i+1)++row(diffLin_right,i+1))
 );

% 线性部分建模
predicate linear_round(
  array [0..1, 0..word_size-1] of var 0..1 : X, 
  array [0..word_size-1] of var 0..1 : p, 
  array [0..1, 0..word_size-1] of var 0..1 : Y,   
  array [0..word_size-1] of var 0..1 : lin_and1_,
  array [0..word_size-1] of var 0..1 : lin_and2_,
  int: round) = (
   let {
     array [0..word_size-1] of var 0..1: X1 = RRot(lin_and1_, 1),  
     array [0..word_size-1] of var 0..1: X2 = RRot(lin_and2_, 8),
     array [0..word_size-1] of var 0..1: X3 = RRot(X[1,0..word_size-1], 2),
   } in
   linear_and_word(lin_and1_,lin_and2_, X[1,0..word_size-1],p) /\%与操作，没有问题
   three_way_fork_word(X[0,0..word_size-1], X1, X2, X3,Y[1,0..word_size-1]) 
 );

constraint forall (i in 1..nr_lin) (
 linear_round(
   lin_state[i-1, 0..1, 0..word_size-1], 
   p_lin[i-1, 0..word_size-1], 
   lin_state[i, 0..1, 0..word_size-1], 
   lin_and1[i-1, 0..word_size-1],
   lin_and2[i-1, 0..word_size-1],
   i-1) /\ lin_state[i, 0, 0..word_size-1] = lin_state[i-1, 1, 0..word_size-1]%注意左右交换
 ) /\ (sum(lin_state[0, 0..1, 0..word_size-1])>0) %输入线性掩码不为0
 /\ (sum(lin_state[nr_lin, 0..1, 0..word_size-1])<=2) /\(lin_state[nr_lin,0,0]==0);


 
% Rotational constants(移位常数)
int: lr1=0; % left rotation value : 1
int: lr8=5; % left rotation value : 8
int: lr2=1; % left rotation value : 2

% (差分部分与差分-线性的链接限制)
constraint forall(i in 0..word_size-1)(
             if left_diff[nr_diff, i] = 0 then 
               diffLin_left[0, i] = -1 
             else 
               diffLin_left[0, i] = 1 
             endif
           );
constraint forall(i in 0..word_size-1)(
             if right_diff[nr_diff, i] = 0 then 
               diffLin_right[0, i] = -1 
             else 
               diffLin_right[0, i] = 1 
             endif
           );

% (线性部分与差分-线性的链接限制)
array [1..2*word_size] of var lower..upper: linear_mask_times_diff_lin_output=array1d(0..word_size-1, [if lin_state[0, 0, i]=0 then 1 else lin_state[0, 0, i]*(abs(row(diffLin_left,nr_diffLin)[i])) endif |i in 0..word_size-1])++array1d(0..word_size-1, [if lin_state[0, 1, i]=0 then 1 else lin_state[0, 1, i]*(abs(row(diffLin_right,nr_diffLin)[i])) endif |i in 0..word_size-1]);%返回值为向量

constraint diffLin_corr = product(linear_mask_times_diff_lin_output);
constraint diffLin_corr != 0.0;%中间部分的相关性不等于0



% equation 6 of the paper(逼近方程)
constraint diffLinComplement_corr =
if diffLin_corr <= 0.001 then
-19931.570*diffLin_corr+29.897
elseif diffLin_corr <= 0.004 /\ diffLin_corr > 0.001 then
-584.962*diffLin_corr+10.135
elseif diffLin_corr <= 0.014 /\ diffLin_corr > 0.004 then
-192.645*diffLin_corr+8.506
elseif diffLin_corr <= 0.053 /\ diffLin_corr > 0.014 then
-50.626*diffLin_corr+6.575
elseif diffLin_corr <= 0.142 /\ diffLin_corr > 0.053 then
-11.87*diffLin_corr+4.483
elseif diffLin_corr <= 0.246 /\ diffLin_corr > 0.142 then
-8.613*diffLin_corr+4.020
elseif diffLin_corr <= 0.595 /\ diffLin_corr > 0.246 then
-3.761*diffLin_corr+2.825
elseif diffLin_corr <= 0.998 /\ diffLin_corr > 0.595 then
-1.444*diffLin_corr+1.444
else
1=1
endif;



% getting the exponent of the probability for the top part
constraint forall (i in 0..nr_diff-1,j in 0..word_size-1)(if (varibits[i,j]=doublebits[i,j]) then CC[i,j]=0 else CC[i,j]=1 endif);
constraint diff_probability = sum(CC[0..nr_diff-1, 0..word_size-1]);

% getting the exponent of the probability for the bottom part
constraint lin_corr = sum(p_lin[0..nr_lin, 0..word_size-1]);

%solve :: float_search(diffLin_left,0.000000000001,smallest,indomain_min,complete)
solve minimize  diffLinComplement_corr+diff_probability+2*lin_corr;

output["Differential Part: Round " ++ show(i) ++ "\n" ++ show(left_diff[i, 0..word_size-1]) ++ show(right_diff[i, 0..word_size -1]) ++ "\n"  | i in 0..nr_diff];
output["Differential-Linear Part: Round " ++ show(i) ++ "\n" ++show(diffLin_left[i, 0..word_size-1]) ++ show(diffLin_right[i, 0..word_size-1]) ++ "\n"  | i in 0..nr_diffLin];
output["Linear Part: Round " ++ show(i) ++ "\n" ++ show(lin_state[i, 0, 0..word_size-1]) ++ show(lin_state[i,1,0..word_size-1]) ++ "\n"  | i in 0..nr_lin];
output["\n"];
output["Differential Probability: " ++ show(diff_probability)];
output["\n"];
output["Differential-Linear Probability: " ++ show(diffLin_corr)];
output["\n"];
output["Linear Probability: " ++ show(lin_corr)];
output["\n"];
output["Total Number of Rounds:" ++ show(nr_diff+nr_diffLin+nr_lin)];
output["\n"];
output["Total Probability:" ++ show(log2(diffLin_corr)-diff_probability-2*lin_corr)];
%output["\n"];
%output["Python Format Distinguisher"];
%output["\n"];
%output["diff_left = int('" ++ int_array_to_bin_str(row(left_diff,0)) ++ "', 2)"];
%output["\n"];
%output["diff_right = int('" ++ int_array_to_bin_str(row(right_diff,0)) ++ "', 2)"];
%output["\n"];
%output["lin_mask_left = int('" ++ int_array_to_bin_str(lin_state[nr_lin, 0, 0..word_size-1]) ++ "', 2)"];
%output["\n"];
%output["lin_mask_right = int('" ++ int_array_to_bin_str(lin_state[nr_lin, 1, 0..word_size-1]) ++ "', 2)"];