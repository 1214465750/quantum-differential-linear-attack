%与操作
predicate linear_and_bit(var 0..1:a,var 0..1:b,var 0..1:c,var 0..1:d) =%d表示相关性的指数的绝对值
 %-d >= -1 /\
 %b >= 0 /\
 %a >= 0 /\
 -b + d >= 0 /\
 -a + d >= 0 /\
 -c + d >= 0 /\
 -d + c >= 0;

predicate linear_and_word(
 array [int] of var 0..1 : A, 
 array [int] of var 0..1 : B, 
 array [int] of var 0..1 : C, 
 array [int] of var 0..1 : di_list) = 
let {
    int: n = length(A);
} in 
 forall (j in 0..n-1) (linear_and_bit(A[j],B[j], C[j],di_list[j]));

%三分支
/*predicate three_way_fork_word(array [int] of var 0..1: X, array [int] of var 0..1: Y, array [int] of var 0..1: Z, array [int] of var 0..1: Dummy) = 
let {
    int: n = length(X);
} in 
forall (j in 0..n-1) (
three_way_fork(X[j], Y[j], Z[j], Dummy[j])
);
 
predicate three_way_fork(var 0..1:x, var 0..1:y, var 0..1:z, var 0..1:dummy) =
  x + y + z - 2*dummy >= 0 /\
  x + y + z <= 2 /\
  dummy - x >= 0 /\
  dummy - y >= 0 /\
  dummy - z >= 0;*/
predicate three_way_fork(var 0..1:x, var 0..1:y, var 0..1:z, var 0..1:w, var 0..1:v) =
  %-y >= -1 /\
  %-x >= -1 /\
  %-z >= -1 /\
  %v >= 0 /\
  %-w >= -1 /\
  %w >= 0 /\
  %z >= 0 /\
  x - y + z - w - v >= -2 /\
  %x >= 0 /\
  x + y - z - w - v >= -2 /\
  %y >= 0 /\
  x - y - z + w - v >= -2 /\
  x + y + z + w - v >= 0 /\
  x + y + z - w + v >= 0 /\
  x + y - z + w + v >= 0 /\
  %-v >= -1 /\
  -x + y + z + w + v >= 0 /\
  x - y + z + w + v >= 0 /\
  -x + y + z - w - v >= -2 /\
  -x + y - z + w - v >= -2 /\
  -x - y + z + w - v >= -2 /\
  -x - y - z - w - v >= -4 /\
  -x - y - z + w + v >= -2 /\
  -x - y + z - w + v >= -2 /\
  x - y - z - w + v >= -2 /\
  -x + y - z - w + v >= -2;
 
predicate three_way_fork_word(array [int] of var 0..1: X, array [int] of var 0..1: Y, array [int] of var 0..1: Z, array [int] of var 0..1: W, array [int] of var 0..1: V) = 
let {
    int: n = length(X);
} in 
forall (j in 0..n-1) (
three_way_fork(X[j], Y[j], Z[j], W[j], V[j])
);
 

 
