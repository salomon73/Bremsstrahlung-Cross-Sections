% compile the hypergeometric2f1_SPECIFYEXTENSION.cpp file

% S. Guinchard EPFL <salomon.guinchard@epfl.ch>

mex -I/opt/homebrew/include ...
    -L/opt/homebrew/lib -lflint -lgmp -lmpfr -lmpc ...
    hypergeometric2f1_flint_vec.cpp