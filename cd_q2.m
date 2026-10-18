% MAIN FUNCTION (note this is also known as the convection-diffusion problem)
%
% num_nodes (Natural number, i.e. integer)
%       Number of nodes after 0 for x and y, e.g. 2 gives x,y={0 , 0.5 , 1}
%           Note this is assuming Q1 Mesh (i.e. for element partitions)

function cd_q2(num_nodes)
    tic;

    % ====================================================================
    %% 0 - Numerical functions for later reference
    % ====================================================================
    % Solutions, external force f, and velocity field beta
    u_x = @(x,y) sin(pi * x) .* sin(pi * y);
    u_y = @(x,y) -sin(pi * x) .* sin(pi * y);
    u_x_dx = @(x,y) pi * cos(pi * x) .* sin(pi * y);
    u_x_dy = @(x,y) pi * sin(pi * x) .* cos(pi * y);
    u_y_dx = @(x,y) -pi * cos(pi * x) .* sin(pi * y);
    u_y_dy = @(x,y) -pi * sin(pi * x) .* cos(pi * y);

    beta_x = @(x,y) 4*(2*y-1).*(1-x).*x;
    beta_y = @(x,y) -4*(2*x-1).*(1-y).*y;

    f_x = @(x,y) 2 * pi^2 * sin(pi * x) .* sin(pi * y) + beta_x(x,y) * u_x_dx(x,y) + beta_y(x,y) * u_x_dy(x,y);
    f_y = @(x,y) -2 * pi^2 * sin(pi * x) .* sin(pi * y) + beta_x(x,y) * u_y_dx(x,y) + beta_y(x,y) * u_y_dy(x,y);


    % Q2 basis functions and gradients, constructed from zero lines
    z1 = @(x) x;
    z2 = @(x) 1-2*x;
    z3 = @(x) 1-x;
    z4 = @(y) y;
    z5 = @(y) 1-2*y;
    z6 = @(y) 1-y;

    zeta1 = @(x,y)       z2(x) .* z3(x) .* z5(y) .* z6(y);
    zeta1i = @(x,y) 4 *  z1(x) .* z3(x) .* z5(y) .* z6(y);
    zeta2 = @(x,y)  -    z1(x) .* z2(x) .* z5(y) .* z6(y);
    zeta2i = @(x,y) -4*  z1(x) .* z2(x) .* z4(y) .* z6(y);
    zeta3 = @(x,y)       z1(x) .* z2(x) .* z4(y) .* z5(y);
    zeta3i = @(x,y) -4 * z1(x) .* z3(x) .* z4(y) .* z5(y);
    zeta4 = @(x,y)  -    z2(x) .* z3(x) .* z4(y) .* z5(y);
    zeta4i = @(x,y) 4 *  z2(x) .* z3(x) .* z4(y) .* z6(y);
    zetam = @(x,y) 16 *  z1(x) .* z3(x) .* z4(y) .* z6(y);

    % Group zetas together under one function for convenience
    Zeta = @(x,y) [ zeta1(x,y) , zeta1i(x,y) , zeta2(x,y) , ...
                        zeta2i(x,y) , zeta3(x,y) , zeta3i(x,y) , ...
                        zeta4(x,y) , zeta4i(x,y) , zetam(x,y)];
    %disp([Zeta(0,0);Zeta(0.5,0);Zeta(1,0);Zeta(1,0.5);Zeta(1,1);Zeta(0.5,1);Zeta(0,1);Zeta(0,0.5);Zeta(0.5,0.5)])

    g_zeta1 = @(x,y)       [z5(y) .* z6(y) .* (-2*z3(x) - z2(x)) ; ...
                            z2(x) .* z3(x) .* (-2*z6(y) - z5(y))];
    g_zeta1i = @(x,y) 4 *  [z5(y) .* z6(y) .* (z3(x) - z1(x)) ; ...
                            z1(x) .* z3(x) .* (-2*z6(y) - z5(y))];
    g_zeta2 = @(x,y)  -    [z5(y) .* z6(y) .* (z2(x) -2* z1(x)) ; ...
                            z1(x) .* z2(x) .* (-2*z6(y) - z5(y))];
    g_zeta2i = @(x,y) -4 * [z4(y) .* z6(y) .* (z2(x) -2* z1(x)) ; ...
                            z1(x) .* z2(x) .* (z6(y) - z4(y))];
    g_zeta3 = @(x,y)       [z4(y) .* z5(y) .* (z2(x) -2* z1(x)) ; ...
                            z1(x) .* z2(x) .* (z5(y) -2* z4(y))];
    g_zeta3i = @(x,y) -4 * [z4(y) .* z5(y) .* (z3(x) - z1(x)) ; ...
                            z1(x) .* z3(x) .* (z5(y) -2* z4(y))];
    g_zeta4 = @(x,y)  -    [z4(y) .* z5(y) .* (-2*z3(x) - z2(x)) ; ...
                            z2(x) .* z3(x) .* (z5(y) -2* z4(y))];
    g_zeta4i = @(x,y) 4 *  [z4(y) .* z6(y) .* (-2*z3(x) - z2(x)) ; ...
                            z2(x) .* z3(x) .* (z6(y) - z4(y))];
    g_zetam = @(x,y)  16 * [z4(y) .* z6(y) .* (z3(x) - z1(x)) ; ...
                            z1(x) .* z3(x) .* (z6(y) - z4(y))];
    % Group gradient zetas together under one function for convenience
    Grad_zeta = @(x,y)[g_zeta1(x,y) , g_zeta1i(x,y) , g_zeta2(x,y), ...
                        g_zeta2i(x,y) , g_zeta3(x,y) , g_zeta3i(x,y), ...
                        g_zeta4(x,y) , g_zeta4i(x,y) , g_zetam(x,y)];


    % Recurring floats saved for optimization
    a = 25/324;
    b = 10/81;
    c = (1-sqrt(3/5)) / 2;
    d = (1+sqrt(3/5)) / 2;
    % Assemble quadrature weights and nodes
    Qweights = [a,b,a,b, 16/81 ,b,a,b,a];
    Qx_nodes = [c,c,c,1/2,1/2,1/2,d,d,d];
    Qy_nodes = [c,1/2,d,c,1/2,d,c,1/2,d];
    % % Recurring floats saved for optimization
    % w1 = ( (18 - sqrt(30))/(36) )^2 / 4;
    % w2 = (49 / 216) / 4;
    % w3 = ( (18 + sqrt(30))/(36) )^2 / 4;
    % x1 = (1 - sqrt(3/7 + 2/7 * sqrt(6/5))) / 2;
    % x2 = (1 - sqrt(3/7 - 2/7 * sqrt(6/5))) / 2;
    % x3 = (1 + sqrt(3/7 - 2/7 * sqrt(6/5))) / 2;
    % x4 = (1 + sqrt(3/7 + 2/7 * sqrt(6/5))) / 2;
    % % Assemble quadrature weights and nodes
    % Qweights = [w1,w2,w2,w1 , w2,w3,w3,w2 , w2,w3,w3,w2 , w1,w2,w2,w1];
    % Qx_nodes = [x1,x1,x1,x1 , x2,x2,x2,x2 , x3,x3,x3,x3 , x4,x4,x4,x4];
    % Qy_nodes = [x1,x2,x3,x4 , x1,x2,x3,x4 , x1,x2,x3,x4 , x1,x2,x3,x4];



   
    % ====================================================================
    %% 2 - Q2 Setup
    % ====================================================================
    % Total number of rectangular elements for *both* Q1 and Q2
    num_elem = num_nodes^2;  
    num_nodes_all2 = 2*num_nodes+1;  % Total nodes per variable
    % List of all distinct values
    xy_vals2 = linspace( 0 , 1 , num_nodes_all2);  
    [X2,Y2] = meshgrid(xy_vals2 , xy_vals2);  % Set up mesh grid
    %disp(X2)
    %disp(Y2)

    % Assign a unique global index to each unique (i,j) in the above meshes
    % Will count up from top to bottom and then from left to right
    total_nodes2 = (num_nodes_all2)^2;  % Total nodes period
    glob_idx2 = reshape( linspace(1, total_nodes2, total_nodes2), ...
        [num_nodes_all2 , num_nodes_all2] );
    %disp(glob_idx2)


    % Create local connectivity array
    % (1,1i,2,2i,3,3i,4,4i,m) --> (g1,g2,g3,g4,g5,g6,g7,g8,g9)

    elem = 0;
    connections2 = zeros(9 , num_elem);
    for i = 1 : 2 : (num_nodes_all2 - 2)
        for j = 1 : 2 : (num_nodes_all2 - 2)
            elem = elem + 1;

            % 00 > 10 > 20 > 21 > 22 > 12 > 02 > 01 > 11
            connections2(:,elem) = [glob_idx2(i , j);  % 1
                                    glob_idx2(i , j+1);  % 1i 
                                    glob_idx2(i , j+2);  % 2
                                    glob_idx2(i+1 , j+2);  % 2i
                                    glob_idx2(i+2 , j+2);  % 3
                                    glob_idx2(i+2 , j+1);  % 3i
                                    glob_idx2(i+2 , j);  % 4
                                    glob_idx2(i+1 , j);  % 4i
                                    glob_idx2(i+1 , j+1)];  % m
        end
    end
    %disp(connections2)
     
  


    % ====================================================================
    %% 3 - Begin the Element Iterations
    % ====================================================================
    % Set up our consistent BK and its results
    h = 1 / num_nodes;
    BK = [h , 0; 0 , h];  % BK itself
    det_BK = h^2; % |det(BK)|
    inv_BK = [(1/h) , 0 ; 0 , (1/h)];  % BK^(-1), same as BK^(-T)


    %We have 2 global stiffness matrices and 2 global load vectors:
    %          A and B; Fx and Fy
    %First we'll set up sparse matrices and zero vectors
    A = sparse(total_nodes2,total_nodes2);
    B = A;
    Fx = zeros(total_nodes2,1);
    Fy = Fx;


    for elem = 1:num_elem
        % Set up local stiffness matrices and load vectors
        lA = zeros(9,9);
        lB = lA;
        lFx = zeros(9,1);
        lFy = lFx;

        % Grab the corresponding 1st local node
        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx2 == connections2(1,elem),1);
        bK = [X2(row,col) ; Y2(row,col)];
        %disp(bK)


        % Set up A and B matrices using Gaussian Quadratures
        for i = 1:9  % Rows
            for j = 1:9  % Columns
                % We'll set up each element of our local stiffness
                % matrices using a Gaussian quadrature
                for q = 1:9  % Quadrature index for weight/node
                    Z = Zeta(Qx_nodes(q) , Qy_nodes(q));
                    Gz = Grad_zeta(Qx_nodes(q) , Qy_nodes(q));
                    Xg = bK + BK*[Qx_nodes(q);Qy_nodes(q)];

                    lA(i,j) = lA(i,j) + Qweights(q) * Gz(:,i)' * inv_BK * inv_BK * Gz(:,j);
                    lB(i,j) = lB(i,j) + Qweights(q) * [beta_x(Xg(1) , Xg(2)) , beta_y(Xg(1) , Xg(2))] * inv_BK * Gz(:,i) * Z(j);
                end
            end
        end
        lA = det_BK * lA;
        lB = det_BK * lB;
        %disp(lA)
        %disp(lB)

        
        % Set up F vectors using Gaussian Quadratures
        for i = 1:9  % Rows
            % We'll set up each element of our local load vectors using
            % a Gaussian quadrature
            for q = 1:9  % Quadrature index for weight/node
                Z = Zeta(Qx_nodes(q) , Qy_nodes(q));
                Xg = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];

                lFx(i) = lFx(i) + Qweights(q) * f_x(Xg(1) , Xg(2)) * Z(i);
                lFy(i) = lFy(i) + Qweights(q) * f_y(Xg(1) , Xg(2)) * Z(i);
            end
        end
        lFx = det_BK * lFx;
        lFy = det_BK * lFy;
        %disp(lFx)
        %disp(lFy)


        % Obtain Q1 and Q2 global indices
        g_ind2 = connections2(: , elem);

        % Update values of global stiffness matrices and load vectors
        A(g_ind2,g_ind2) = A(g_ind2,g_ind2) + lA;
        B(g_ind2,g_ind2) = B(g_ind2,g_ind2) + lB;
        Fx(g_ind2) = Fx(g_ind2) + lFx;
        Fy(g_ind2) = Fy(g_ind2) + lFy;
    end




    % ====================================================================
    %% 4 - Set up system, implement boundary conditions, and solve
    % ====================================================================
    % Combined stiffness matrix K will be of the form:
    K = [A+B' , zeros(total_nodes2,total_nodes2);
         zeros(total_nodes2,total_nodes2), A+B'];
    % Meanwhile the combined load vector g will be:
    g = [Fx; Fy];
    %disp(full(K))
    %disp(g)


    % Create 1D array of all Q2 global indices on the boundary
    bound_idx = [];
    bound_idx = union(bound_idx , glob_idx2( 1 , : ));
  	bound_idx = union(bound_idx , glob_idx2( num_nodes_all2 , : ));
    bound_idx = union(bound_idx , glob_idx2( : , 1 ));
  	bound_idx = union(bound_idx , glob_idx2( : , num_nodes_all2 ));
    % HOWEVER, using the above indices alone would only count for u_x's
    % boundary conditions, so to also factor in u_y's conditions, we
    % would need to the same indices translated the length of A+Bxx,
    % which happens to be num_nodes_all2.
    bound_idx = union(bound_idx, bound_idx + total_nodes2);
    %disp(bound_idx)


    % Set all rows and columns corresponding to bound indices to 0
    %spy(K)
    %rank(full(K))
    K(bound_idx , :) = 0;
    K(: , bound_idx) = 0;
    g(bound_idx) = 0;
    % Update corresponding diagonals of K to 1
    for idx = bound_idx(:)'
        K(idx , idx) = 1;
    end


    %disp(full(K))
    %spy(K)
    %disp(rank(full(K)))
    %disp(size(full(K)))
    %disp(g)


    % Solve the system Kz=g - Note that much faster solvers are viable,
    % but we'll save those for once this script works.
    z = K \ g;
    %disp(z)
    uh_x = z(1:total_nodes2);
    uh_y = z(total_nodes2+1:2*total_nodes2);




    % ====================================================================
    %% 5 - Output results
    % ====================================================================
    % Create arrays of our exact solutions
    % Ux = [];
    % Uy = [];
    % for i = 1:length(xy_vals2)
    %     Ux = [Ux , u_x(xy_vals2(i),xy_vals2)];
    %     Uy = [Uy , u_y(xy_vals2(i),xy_vals2)];
    % end
    % Ux = Ux';
    % Uy = Uy';
    % 
    % % Generate and print L2 error norms
    % l2_norm_ux = norm(abs(Ux - uh_x)) / (total_nodes2-1);
    % l2_norm_uy = norm(abs(Uy - uh_y)) / (total_nodes2-1);
    % fprintf('Discrete ux L2 Norm = %12.12f\n', l2_norm_ux);
    % fprintf('Discrete uy L2 Norm = %12.12f\n', l2_norm_uy);


    %toc  % Output runtime
    %fprintf('\n') 



    
    % ====================================================================
    %% 6 - (NEW) Generate and print more accurate L2 and H1 error estimates
    % ====================================================================
    % This'll be the result of summing the L2/H1 errors over the elements,
    % so we first set the total accumulated errors to 0
    ux_l2 = 0;
    uy_l2 = 0;
    ux_h1 = 0;
    uy_h1 = 0;


    % Now we loop over each element to obtain individual L2 and H1 errors
    for elem = 1:num_elem
        n1 = connections2(1,elem);
        n1i = connections2(2,elem);
        n2 = connections2(3,elem);
        n2i = connections2(4,elem);
        n3 = connections2(5,elem);
        n3i = connections2(6,elem);
        n4 = connections2(7,elem);
        n4i = connections2(8,elem);
        nm = connections2(9,elem);

        % Grab the corresponding 1st local node (same for both Q2 and Q1)
        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx2 == n1,1);
        bK = [X2(row,col) ; Y2(row,col)];

        % On the reference element:
        %   uh = uh(n1) * zeta1 + uh(n1i) * zeta1i + ...
        UH_x = @(x,y) uh_x(n1) * zeta1(x,y) + uh_x(n1i) * zeta1i(x,y) + uh_x(n2) * zeta2(x,y) + ...
            uh_x(n2i) * zeta2i(x,y) + uh_x(n3) * zeta3(x,y) + uh_x(n3i) * zeta3i(x,y) + ...
            uh_x(n4) * zeta4(x,y) + uh_x(n4i) * zeta4i(x,y) + uh_x(nm) * zetam(x,y);
        UH_y = @(x,y) uh_y(n1) * zeta1(x,y) + uh_y(n1i) * zeta1i(x,y) + uh_y(n2) * zeta2(x,y) + ...
            uh_y(n2i) * zeta2i(x,y) + uh_y(n3) * zeta3(x,y) + uh_y(n3i) * zeta3i(x,y) + ...
            uh_y(n4) * zeta4(x,y) + uh_y(n4i) * zeta4i(x,y) + uh_y(nm) * zetam(x,y);
        % Meanwhile on the reference element:
        %   gradient uh = uh(n1) * g_zeta1 + uh(n1i) * g_zeta1i + ...
        g_UH_x = @(x,y) uh_x(n1) * g_zeta1(x,y) + uh_x(n1i) * g_zeta1i(x,y) + uh_x(n2) * g_zeta2(x,y) + ...
            uh_x(n2i) * g_zeta2i(x,y) + uh_x(n3) * g_zeta3(x,y) + uh_x(n3i) * g_zeta3i(x,y) + ...
            uh_x(n4) * g_zeta4(x,y) + uh_x(n4i) * g_zeta4i(x,y) + uh_x(nm) * g_zetam(x,y);
        g_UH_y = @(x,y) uh_y(n1) * g_zeta1(x,y) + uh_y(n1i) * g_zeta1i(x,y) + uh_y(n2) * g_zeta2(x,y) + ...
            uh_y(n2i) * g_zeta2i(x,y) + uh_y(n3) * g_zeta3(x,y) + uh_y(n3i) * g_zeta3i(x,y) + ...
            uh_y(n4) * g_zeta4(x,y) + uh_y(n4i) * g_zeta4i(x,y) + uh_y(nm) * g_zetam(x,y);
        % Note that we should also multiply inv_BK to both gradient
        % functions to reflect the transformation from K to K-hat, but
        % we'll do this in the upcoming lines where we update the H1 norms


        for q = 1:9  % Quadrature index for weight/node
            xr = Qx_nodes(q);  % x in reference element
            yr = Qy_nodes(q);  % y in reference element
            G_node = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];
            xg = G_node(1);  % Global x in current element
            yg = G_node(2);  % Global y in current element

            % Update L2 norms of ux and uy
            ux_l2 = ux_l2 + det_BK * Qweights(q) * (UH_x(xr,yr) - u_x(xg,yg))^2;
            uy_l2 = uy_l2 + det_BK * Qweights(q) * (UH_y(xr,yr) - u_y(xg,yg))^2;

            % Update H1 norms of ux and uy
            ux_h1 = ux_h1 + det_BK * Qweights(q) * norm(inv_BK * g_UH_x(xr,yr) - [u_x_dx(xg,yg) ; u_x_dy(xg,yg)])^2;
            uy_h1 = uy_h1 + det_BK * Qweights(q) * norm(inv_BK * g_UH_y(xr,yr) - [u_y_dx(xg,yg) ; u_y_dy(xg,yg)])^2;
        end
    end


    % Update our errors since they're the square roots of the above sums
    ux_l2 = sqrt(ux_l2);
    uy_l2 = sqrt(uy_l2);
    ux_h1 = sqrt(ux_h1);
    uy_h1 = sqrt(uy_h1);

    % Print the results
    fprintf('ux L2 Norm = %12.12f\n', ux_l2);
    fprintf('uy L2 Norm = %12.12f\n', uy_l2);
    %fprintf('ux H1 Norm = %12.12f\n', ux_h1);
    %fprintf('uy H1 Norm = %12.12f\n', uy_h1);
    %fprintf('\n') 
end