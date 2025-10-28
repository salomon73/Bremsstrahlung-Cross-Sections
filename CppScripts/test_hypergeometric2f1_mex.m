function test_hypergeometric2f1_mex()
    % Test the MEX wrapper for FLINT hypergeometric 2F1 (scalar)

    fprintf('Running unit test for hypergeometric2f1_mex...\n');

    % Test cases: scalar inputs
    test_cases = {
        {0.5, 1.0, 1.5, 0.25},
        {1.0, 2.0, 3.0, 0.5},
        {0.2, 0.3, 0.5, -0.7}
    };

    % Compare with MATLAB built-in hypergeom function (only for real numbers)
    for k = 1:length(test_cases)
        a = test_cases{k}{1};
        b = test_cases{k}{2};
        c = test_cases{k}{3};
        z = test_cases{k}{4};

        % Call MEX
        y_mex = hypergeometric2f1_mex(a, b, c, z);

        % Call MATLAB built-in hypergeom
        y_matlab = hypergeom([a,b], c, z);

        % Display results
        fprintf('Test %d: 2F1(%g,%g;%g;%g)\n', k, a, b, c, z);
        fprintf('   MEX result   : %g + %gi\n', real(y_mex), imag(y_mex));
        fprintf('   MATLAB result: %g + %gi\n', real(y_matlab), imag(y_matlab));
        fprintf('   Relative error: %.2e\n', abs(y_mex - y_matlab)/abs(y_matlab));
        fprintf('\n');
    end

    fprintf('Unit test finished.\n');
end
