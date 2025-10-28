function test_hypergeometric2f1_mex_vec(n,m)
    disp('Running unit test for hypergeometric2f1_flint_vec...');

    % Test input arrays
    a = 0.5*ones(1,m);
    b = 1*ones(1,m);
    c = 1.5*ones(1,m);
    z = 0.25*ones(1,m);

    % Call vectorized MEX
    tic,
    y_mex = hypergeometric2f1_flint_vec(a, b, c, z);
    t1 = toc;
    if any(y_mex-y_mex(1)), error('hypergeometric2f1_flint_vec not working as intended'), end
    
    % Compute MATLAB 
    tic,
    y_matlab = arrayfun(@(aa,bb,cc,zz) hypergeom([aa,bb],cc,zz), a,b,c,z);  % built-in function
    t2 = toc;
    if any(y_matlab-y_matlab(1)), error('hypergeom not working as intended'), end

    % Compute MATLAB hypergeometric2f1
    tic,
    y_matlab2f1 = hypergeometric2f1(a,b,c,z,n);
    t3 = toc;
    if any(y_matlab2f1-y_matlab2f1(1)), error('hypergeometric2f1 not working as intended'), end

    % Display results
    fprintf('Test 2F1(%g,%g;%g;%g)\n', a(1), b(1), c(1), z(1));
    fprintf('   MEX    : %.2f s; res: %.6f + %.6fi\n', t1, real(y_mex(1)), imag(y_mex(1)));
    fprintf('   MATLAB : %.2f s; res: %.6f + %.6fi\n', t2, real(y_matlab(1)), imag(y_matlab(1)));
    fprintf('   MAT2F1 : %.2f s; res: %.6f + %.6fi\n', t3, real(y_matlab2f1(1)), imag(y_matlab2f1(1)));
    relerr = abs(y_mex(1) - y_matlab(1)) / max(abs(y_matlab(1)), eps);
    relerr2 = abs(y_matlab2f1(1) - y_matlab(1)) / max(abs(y_matlab(1)), eps);
    fprintf('   Relative error mex/mat: %.2e\n\n', relerr);
    fprintf('   Relative error 2f1/mat: %.2e\n\n', relerr2);

    disp('Unit test finished.');
end
