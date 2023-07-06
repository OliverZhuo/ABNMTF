function [U,W,out] = ABNMTF(A, k)
epsilon = sqrt(eps);

n = size(A,1);
U = rand(n,k);
W = eye(k,k);
O = eye(k);
B = eye(k);

a = 1;
b = 100;
c = 1;
d = 100;
lambda = 100;

n_iter = 500;


for iter = 1:n_iter
    U = U .* ((A'+A)*U*(W+W')) ./ (2.*(U*W'*U'*U*W'+U*W*U'*U*W+U*O ) + epsilon);
    W = W .* (U'*A*U+U'*A'*U+2*lambda.*eye(k,k) ) ./ (2.*U'*U*W*U'*U+W*B+2*lambda.*W + epsilon);

    alpha = (n + a - 1) ./ (sum(U.*U,1)+b);
    beta = (0.5*n + c - 1) ./ (0.5*sum(W.*W,1)+d);
    for i=1:k
        O(i,i) = alpha(i);
        B(i,i) = beta(i);
    end

end

end
