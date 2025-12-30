float:lower=-1;
float:upper=1;
%异或

function var lower..upper: continuous_xor_bit(var lower..upper: x, var lower..upper: y) = -1*x*y;


function array [int] of var lower..upper: continuous_xor(array [int] of var float: aa,array [int] of var float: bb) = 
let {
    int: n = length(aa);
} in 
array1d(0..n-1, [continuous_xor_bit(aa[i], bb[i])|i in 0..n-1]);

%与操作
function var lower..upper: continuous_and_bit(var lower..upper: x, var lower..upper: y) = 0.25*(x+y-x*y-1);

function array [int] of var lower..upper: continuous_and(array [int] of var float: aa,array [int] of var float: bb) = 
let {
    int: n = length(aa);
} in 
array1d(0..n-1, [continuous_and_bit(aa[i], bb[i])|i in 0..n-1]);


%循环左移
function array[int] of var lower..upper: continuous_LRot(array[int] of var lower..upper: X, int: val)=
let {
    int: n = length(X);
} in 
array1d(0..n-1, [X[(j+val) mod n] | j in 0..n-1]);

%差分部分与差分-线性部分的连接限制
function array [int] of var lower..upper: cast(array [int] of var 0..1:int_var) = 
let {
    int: n = length(int_var);
} in 
    array1d(0..n-1, [if int_var[i] = 0 then -1.0 else 1.0 endif|i in 0..n-1]);