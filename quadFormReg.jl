function denseForm(A, B, N, Q, R, umin, umax, xmin, xmax, x0)
	
	# Model~
	nx, nu = size(B);

	At = eye(nx,nx);
    for i = 1:(N-1)
    	At = [At; A^i];
    end

    Bt = zeros(N*nx,(N-1)*nu);
	bt = B;
	for j = N-2:-1:1
		bt = [A^(N-2-j+1)*B bt]
	end
	for i=nx:nx:((N-1)*nx)
		Bt[i+1:i+nx,1:nu*round(Int,(i/nx))] = bt[:,round(Int,(N-1)*nu-(i/nx*nu)+1):end];
	end

	# Q~, R~
	Qt = kron(eye(N), Q);
	Rt = kron(eye(N), R);

	# constraints~
	Umin = umin;
	Umax = umax;
	Xmin = xmin;
	Xmax = xmax;

	for i = 1:N-1
		Umin = [Umin; umin];
		Umax = [Umax; umax];
		Xmin = [Xmin; xmin];
		Xmax = [Xmax; xmax];
	end

	G = [-eye(N*nu,(N-1)*nu), eye(N*nu,(N-1)*nu), -Bt, Bt];
	H =  [-Umin, Umax, -Xmin + At*x0; Xmax - At*x0];

	P = 2*Bt'*Qt*Bt + Rt;
	q = 2*x0'*At'*Qt*Bt;
	r = x0'*At'*Qt*At*x0;

	return P, q, r, G, H
end