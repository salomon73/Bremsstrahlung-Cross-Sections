function test_hypergeometric2f1_mex_vec()
    disp('Running unit test for hypergeometric2f1_flint_vec...');

    % Test input arrays
    a = [0.5, 1, 0.2];
    b = [1, 2, 0.3];
    c = [1.5, 3, 0.5];
    z = [0.25, 0.5, -0.7];

    % Call vectorized MEX
    y_mex = hypergeometric2f1_flint_vec(a, b, c, z);

    % Compute MATLAB results
    y_matlab = arrayfun(@(aa,bb,cc,zz) hypergeom([aa,bb],cc,zz), a,b,c,z);  % built-in function

    % Reshape for comparison
    y_matlab = reshape(y_matlab, size(a));

    % Display results
    for i = 1:length(a)
        fprintf('Test %d: 2F1(%g,%g;%g;%g)\n', i, a(i), b(i), c(i), z(i));
        fprintf('   MEX result   : %.6f + %.6fi\n', real(y_mex(i)), imag(y_mex(i)));
        fprintf('   MATLAB result: %.6f + %.6fi\n', real(y_matlab(i)), imag(y_matlab(i)));
        relerr = abs(y_mex(i) - y_matlab(i)) / max(abs(y_matlab(i)), eps);
        fprintf('   Relative error: %.2e\n\n', relerr);
    end

    disp('Unit test finished.');
end
