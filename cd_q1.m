% MAIN FUNCTION (note this is also known as the convection-diffusion problem)
%
% num_nodes (Natural number, i.e. integer)
%       Number of nodes after 0 for x and y, e.g. 2 gives x,y={0 , 0.5 , 1}
%           Note this is assuming Q1 Mesh (i.e. for element partitions)

function cd_q1(num_nodes)
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


     % Q1 basis functions
    eta1 = @(x,y) (1-x) .* (1-y);
    eta2 = @(x,y) x .* (1-y);
    eta3 = @(x,y) x .* y;
    eta4 = @(x,y) (1-x) .* y;
    % Group etas together under one function for convenience
    Eta = @(x,y) [ eta1(x,y) , eta2(x,y) , eta3(x,y) , eta4(x,y)];


    g_eta1 = @(x,y) [y-1 ; x-1];
    g_eta2 = @(x,y) [1-y ; -x];
    g_eta3 = @(x,y) [y ; x];
    g_eta4 = @(x,y) [-y ; 1-x];
    % Group gradient etas together under one function for convenience
    Grad_eta = @(x,y) [g_eta1(x,y) , g_eta2(x,y), ...
                        g_eta3(x,y) , g_eta4(x,y)];


    % Recurring floats saved for optimization
    a = 25/324;
    b = 10/81;
    c = (1-sqrt(3/5)) / 2;
    d = (1+sqrt(3/5)) / 2;
    % Assemble quadrature weights and nodes
    Qweights = [a,b,a,b, 16/81 ,b,a,b,a];
    Qx_nodes = [c,c,c,1/2,1/2,1/2,d,d,d];
    Qy_nodes = [c,1/2,d,c,1/2,d,c,1/2,d];



   
    % ====================================================================
    %% 1 - Q1 Setup
    % ====================================================================
    % Total number of rectangular elements for Q1
    num_elem = num_nodes^2;  
    num_nodes_all1 = num_nodes+1;  % Total nodes per variable in Q1
    %List of all distinct values
    xy_vals1 = linspace( 0 , 1 , num_nodes_all1);  
    [X1,Y1] = meshgrid(xy_vals1 , xy_vals1);  % Set up mesh grid
    %disp(X1)
    %disp(Y1)

    % Assign a unique global index to each unique (i,j) in the above meshes
    % Will count up from top to bottom and then from left to right
    total_nodes1 = (num_nodes_all1)^2;  % Total nodes period
    glob_idx1 = reshape( linspace(1, total_nodes1, total_nodes1), ...
        [num_nodes_all1 , num_nodes_all1] );
    %disp(glob_idx1)


    % Create local connectivity array
    % (1,2,3,4) --> (g1,g2,g3,g4)
    elem = 0;
    connections1 = zeros(4 , num_elem);
    for i = 1:num_nodes
        for j = 1:num_nodes
            elem = elem + 1;

            % Starting at (xi,yj), 00 > 10 > 11 > 01
            connections1(:,elem) = [glob_idx1(i , j);  % 1
                glob_idx1(i , j+1);  % 2
                glob_idx1(i+1 , j+1);  % 3
                glob_idx1(i+1 , j)];  % 4
        end
    end
    %disp(connections1)
     
  


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
    A = sparse(total_nodes1,total_nodes1);
    B = A;
    Fx = zeros(total_nodes1,1);
    Fy = Fx;


    for elem = 1:num_elem
        % Set up local stiffness matrices and load vectors
        lA = zeros(4,4);
        lB = lA;
        lFx = zeros(4,1);
        lFy = lFx;

        % Grab the corresponding 1st local node
        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx1 == connections1(1,elem),1);
        bK = [X1(row,col) ; Y1(row,col)];
        %disp(bK)


        % Set up A and B matrices using Gaussian Quadratures
        for i = 1:4  % Rows
            for j = 1:4  % Columns
                % We'll set up each element of our local stiffness
                % matrices using a Gaussian quadrature
                for q = 1:9  % Quadrature index for weight/node
                    E = Eta(Qx_nodes(q) , Qy_nodes(q));
                    Ge = Grad_eta(Qx_nodes(q) , Qy_nodes(q));
                    Xg = bK + BK*[Qx_nodes(q);Qy_nodes(q)];

                    lA(i,j) = lA(i,j) + Qweights(q) * Ge(:,i)' * inv_BK * inv_BK * Ge(:,j);
                    lB(i,j) = lB(i,j) + Qweights(q) * [beta_x(Xg(1) , Xg(2)) , beta_y(Xg(1) , Xg(2))] * inv_BK * Ge(:,i) * E(j);
                end
            end
        end
        lA = det_BK * lA;
        lB = det_BK * lB;
        %disp(lA)
        %disp(lB)

        
        % Set up F vectors using Gaussian Quadratures
        for i = 1:4  % Rows
            % We'll set up each element of our local load vectors using
            % a Gaussian quadrature
            for q = 1:9  % Quadrature index for weight/node
                E = Eta(Qx_nodes(q) , Qy_nodes(q));
                Xg = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];

                lFx(i) = lFx(i) + Qweights(q) * f_x(Xg(1) , Xg(2)) * E(i);
                lFy(i) = lFy(i) + Qweights(q) * f_y(Xg(1) , Xg(2)) * E(i);
            end
        end
        lFx = det_BK * lFx;
        lFy = det_BK * lFy;
        %disp(lFx)
        %disp(lFy)


        % Obtain Q1 global indices
        g_ind1 = connections1(: , elem);

        % Update values of global stiffness matrices and load vectors
        A(g_ind1,g_ind1) = A(g_ind1,g_ind1) + lA;
        B(g_ind1,g_ind1) = B(g_ind1,g_ind1) + lB;
        Fx(g_ind1) = Fx(g_ind1) + lFx;
        Fy(g_ind1) = Fy(g_ind1) + lFy;
    end




    % ====================================================================
    %% 4 - Set up system, implement boundary conditions, and solve
    % ====================================================================
    % Combined stiffness matrix K will be of the form:
    K = [A+B' , zeros(total_nodes1,total_nodes1);
         zeros(total_nodes1,total_nodes1), A+B'];
    % Meanwhile the combined load vector g will be:
    g = [Fx; Fy];
    %disp(full(K))
    %disp(g)


    % Create 1D array of all Q1 global indices on the boundary
    bound_idx = [];
    bound_idx = union(bound_idx , glob_idx1( 1 , : ));
  	bound_idx = union(bound_idx , glob_idx1( num_nodes_all1 , : ));
    bound_idx = union(bound_idx , glob_idx1( : , 1 ));
  	bound_idx = union(bound_idx , glob_idx1( : , num_nodes_all1 ));
    % HOWEVER, using the above indices alone would only count for u_x's
    % boundary conditions, so to also factor in u_y's conditions, we
    % would need to the same indices translated the length of A+Bxx,
    % which happens to be num_nodes_all1.
    bound_idx = union(bound_idx, bound_idx + total_nodes1);
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
    uh_x = z(1:total_nodes1);
    uh_y = z(total_nodes1+1:2*total_nodes1);




    % ====================================================================
    %% 5 - Output results
    % ====================================================================
    % Create arrays of our exact solutions
    % Ux = [];
    % Uy = [];
    % for i = 1:length(xy_vals1)
    %     Ux = [Ux , u_x(xy_vals1(i),xy_vals1)];
    %     Uy = [Uy , u_y(xy_vals1(i),xy_vals1)];
    % end
    % Ux = Ux';
    % Uy = Uy';
    % 
    % % Generate and print L2 error norms
    % l2_norm_ux = norm(abs(Ux - uh_x)) / (total_nodes1-1);
    % l2_norm_uy = norm(abs(Uy - uh_y)) / (total_nodes1-1);
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
        n1 = connections1(1,elem);
        n2 = connections1(2,elem);
        n3 = connections1(3,elem);
        n4 = connections1(4,elem);

        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx1 == n1,1);
        bK = [X1(row,col) ; Y1(row,col)];

        UH_x = @(x,y) uh_x(n1) * eta1(x,y) + uh_x(n2) * eta2(x,y) + ...
            uh_x(n3) * eta3(x,y) + uh_x(n4) * eta4(x,y);
        UH_y = @(x,y) uh_y(n1) * eta1(x,y) + uh_y(n2) * eta2(x,y) + ...
            uh_y(n3) * eta3(x,y) + uh_y(n4) * eta4(x,y);
        g_UH_x = @(x,y) uh_x(n1) * g_eta1(x,y) + uh_x(n2) * g_eta2(x,y) + ...
            uh_x(n3) * g_eta3(x,y) + uh_x(n4) * g_eta4(x,y);
        g_UH_y = @(x,y) uh_y(n1) * g_eta1(x,y) + uh_y(n2) * g_eta2(x,y) + ...
            uh_y(n3) * g_eta3(x,y) + uh_y(n4) * g_eta4(x,y);

        for q = 1:9  % Quadrature index for weight/node
            xr = Qx_nodes(q);  % x in reference element
            yr = Qy_nodes(q);  % y in reference element
            G_node = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];
            xg = G_node(1);  % Global x in current element
            yg = G_node(2);  % Global y in current element

            ux_l2 = ux_l2 + det_BK * Qweights(q) * (UH_x(xr,yr) - u_x(xg,yg))^2;
            uy_l2 = ux_l2 + det_BK * Qweights(q) * (UH_y(xr,yr) - u_y(xg,yg))^2;
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