function x = softthre(a, tau)

% solve the optimal problem as: 
%                   argmin 1/2||X - A||_2 + tau||X||_1


x = sign(a).* max( abs(a) - tau, 0);
end