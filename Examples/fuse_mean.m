function [xfkk, X] = fuse_mean(y_a, y_b, X, X_a, X_b, Xij, xf, xab_size, T1, T2, ia1, ia2, cx);
% [xfkk, X, xf] = fuse_mean(y_a, y_b, X, X_a, X_b, Xij, xf, xab_size, T1, T2, ia1, ia2, cx);
% xlkk      Local State vector to be fused
% Slkk      Local Co-Variance Matrix to be fused
% x         Local State vector index
% hr        Number of Local Nodes

x_tmp       = xab_size;

x_a         = T1*y_a;
x_b         = T2*y_b;

x_a(ia1)    = x_b(ia1);
x_b(ia2)    = x_a(ia2);

xij     = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij) + 2*cx*eye(size(x_tmp)))...
    *((pinv(X_b) - pinv(Xij) + cx*eye(size(x_tmp)))*x_a + (pinv(X_a) - pinv(Xij) + cx*eye(size(x_tmp)))*x_b);

xfkk    = X*(pinv(X_a)*x_a + pinv(X_b)*x_b - pinv(Xij)*xij);