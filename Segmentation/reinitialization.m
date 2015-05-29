function phi = reinitialization(phi, dt)

for i=1:10
    old = phi;
    temp = padarray(phi,[1,1],'symmetric','post');
    temp = padarray(temp,[1,1],'symmetric','pre');

    a = phi - temp(1:end-2 , 2:end-1);
    b = temp(3:end , 2:end-1) - phi;
    c = phi - temp(2:end-1 , 1:end-2);
    d = temp(2:end-1 , 3:end) - phi;

    a_p = max(a,0);
    a_m = min(a,0);
    b_p = max(b,0);
    b_m = min(b,0);
    c_p = max(c,0);
    c_m = min(c,0);
    d_p = max(d,0);
    d_m = min(d,0);

    G = zeros(size(phi));
    pos = find(phi>0);
    neg = find(phi<0);
    G(pos) = sqrt(max(a_p(pos).^2 , b_m(pos).^2) + max(c_p(pos).^2 , d_m(pos).^2)) - 1;
    G(neg) = sqrt(max(a_m(neg).^2 , b_p(neg).^2) + max(c_m(neg).^2 , d_p(neg).^2)) - 1;
    sign = phi./abs(phi + eps);
    phi = phi - dt*sign.*G;

    % Stop iteration
    ind = find(abs(phi)<=1);
    M = length(ind);
    Q = sum(abs(phi(ind) - old(ind)))./M;
    if Q <= dt*1^2
        return;
    end
end