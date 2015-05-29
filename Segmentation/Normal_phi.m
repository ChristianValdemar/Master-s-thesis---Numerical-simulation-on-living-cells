function grad_phi = Normal_phi(phi)
% Returns the unit normal vector to phi

[m,n]=size(phi);
%grad_phi = [1:n,1:m,1:2];

phix = [phi(1:m,2)-phi(1:m,1) 1/2*(phi(1:m,3:n)-phi(1:m,1:n-2)) phi(1:m,n)-phi(1:m,n-1)];
phiy = [phi(2,1:n)-phi(1,1:n);1/2*(phi(3:m,1:n)-phi(1:m-2,1:n));phi(m,1:n)-phi(m-1,1:n)];

length_grad_phi = sqrt(phix.^2 + phiy.^2);

grad_phi(:,:,1) = phix./length_grad_phi;
grad_phi(:,:,2) = phiy./length_grad_phi;

