A = [
  3    9   6;
  4   12  12;
  1   -1   1;
];

b = [ 12; 12; 1 ];

[L U P] = lu(A);
y = L \ (P * b);
x = U \ y;
invA = inv(A);
# (ALTERNATIVE) x = invA * b

L
U
P
y
x
invA
