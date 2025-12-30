%循环与操作

predicate SIMON_like_round_f(array [int] of var 0..1 : A, array [int] of var 0..1 : varibits,array [int] of var 0..1 : doublebits, array [int] of var 0..1 : Belt, array [int] of var 0..1 : Z) =
let {
    int: n = length(A);
    int: a = 5;
    int: b = 0;
    int: c = 1;
    int: t1 = 10;
    int: t2 = 5;
} in 
forall (i in 0..n-1) (
-varibits[i] + A[(i+a) mod n] + A[(i+b) mod n] >= 0 /\
varibits[i] - A[(i+b) mod n] >= 0 /\
varibits[i] - A[(i+a) mod n] >= 0 /\%第一个约束
-doublebits[i] - A[(i+a) mod n] >= -1 /\
-doublebits[i] + A[(i+t1) mod n] >= 0 /\
-doublebits[i] + A[(i+b) mod n] >= 0 /\
doublebits[i] - A[(i+b) mod n] + A[(i+a) mod n] - A[(i+t1) mod n] >= -1 /\%第二个约束
A[(i+c) mod n] + Belt[i] - Z[i] >= 0 /\
A[(i+c) mod n] - Belt[i] + Z[i] >= 0 /\
-A[(i+c) mod n] + Belt[i] + Z[i] >= 0 /\
-A[(i+c) mod n] - Belt[i] - Z[i] >= -2 /\%第三个约束
-Z[i] + varibits[i] >= 0 /\ %第四个约束
Z[i] - Z[(i+t2) mod n] - doublebits[i] >= -1 /\
-Z[i] + Z[(i+t2) mod n] - doublebits[i] >= -1);

predicate Xor_bit(var 0..1:a,var 0..1:b,var 0..1:c,var 0..1:d) =
d >= a /\
d >= b /\
d >= c /\
a + b + c >= 2*d /\
a + b + c <= 2;

predicate Xor_word(array [int] of var 0..1 : A, array [int] of var 0..1 : B, array [int] of var 0..1 : C, array [int] of var 0..1 : D, int: r) =
let {
    int: n = length(A);
} in 
forall (j in 0..n-1) (
 Xor_bit(A[j],B[j],C[j], D[j]));
 

  

%function var 0..1: eq(var 0..1:a,var 0..1:b,var int:c) = (if (a=b /\ b=c) then 0 else 1 endif);



 
% Left rotation of X by val positions
function array[int] of var 0..1: LRot(array[int] of var 0..1: X, int: val)=
let {
    int: n = length(X);
} in 
array1d(0..n-1, [X[(j+val) mod n] | j in 0..n-1]);
function array[int] of var 0..1: RRot(array[int] of var 0..1: X, int: val)=
let {
    int: n = length(X);
} in 
array1d(0..n-1, [X[(n+j-val) mod n] | j in 0..n-1]);