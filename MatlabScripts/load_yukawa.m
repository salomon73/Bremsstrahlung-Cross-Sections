function YukawaModel = load_yukawa(dir,Zs, nion, YKorder, plotfig)
%
% This function generates the Multi-Yukawa (MY) model parameters for 
% the screening of the BH (1st Born approximation) bremsstrahlung 
% differential cross-section. 
% 
% see Y. Savoye-Peysson et al 2023 Nucl. Fusion 63 126041 for model.
%
%   Input: 
%
%       - dir     : directory containing the Yukawa DHFS coefficients 
%       - Zs      : Atomic number of the element
%       - nion    : slice of ionization levels to compute the form factor
%       - YKorder : order of the MY model (number of exponentials)
%       - plotfig : bool specifying wether the AFF has to be plotted
%
%   Output:
%   
%       - struct YukawaModel with fields
%
%         - A      : weights of Yukawa exponentials
%         - Lambda : inverse of characteristic length
%         - AFF    : atomic form factor (for postprocessing plot)
%         - qbar   : normalized recoil momentum (for postprocessing plot)
%
%
%
%   Dependencies and toolboxes:
%       
%      Dependencies: 
%
%       - MYi_Zj.m  % DFT multi-Yukawa coefficients        
%           
%   Toolboxes:
%
%       - MATLAB
%
%
%
% Author: Salomon Guinchard 
% <salomon.guinchard@epfl.com>
% Last update: Jul. 17th, 2025.
%
    
    % Add directory path 
    addpath(genpath(dir))
    alpha = 1/137;             % (fine structure constant)
    qbar = logspace(-4,2,400); % Recoil momentum 
    

    if nargin < 2 
	    error(2,'Wrong number of input arguments for load_yukawa');
    end
    if nargin < 3
       nion = [0];
    end
    %
    if nargin < 4
        YKorder = 1;
    end
    if nargin <5
        plotfig = false;
    end

    if Zs == 79
        if YKorder ==4 
            M = MY4_Gold();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_Gold();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_Gold();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_Gold();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end
    elseif Zs==13
        if YKorder ==4 
            M = MY4_Al();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_Al();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_Al();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_Al();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end

        elseif Zs==74
        if YKorder ==4 
            M = MY4_W();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_W();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_W();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_W();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end
        elseif Zs==50
        if YKorder ==4 
            M = MY4_Sn();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_Sn();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_Sn();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_Sn();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end
        elseif Zs==29
        if YKorder ==4 
            M = MY4_Cu();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_Cu();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_Cu();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_Cu();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end
        elseif Zs==4
        if YKorder ==4 
            M = MY4_Be();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:7);
            YukawaModel.Lambda  = neut_coeff(:,8:11);
    
        elseif YKorder ==3
            M = MY3_Be();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:6);
            YukawaModel.Lambda  = neut_coeff(:,7:9);
        
        elseif YKorder ==2
            M = MY2_Be();
            neut_coeff = M(:,:);
            YukawaModel.A = neut_coeff(:,4:5);
            YukawaModel.Lambda  = neut_coeff(:,6:7);
    
        elseif YKorder ==1
            M = MY1_Be();
            neut_coeff = M(:,:);
            YukawaModel.A  = neut_coeff(:,3);
            YukawaModel.Lambda  = neut_coeff(:,4);
    
        else 
          error('YKorder must be 1, 2, 3, or 4.');  
        end
    end

    Ffactor = zeros(numel(nion),length(qbar)); % Atomic Form Factor 
    for i = 1:numel(nion)
        A = YukawaModel.A(nion(i)+1,1:YKorder);
        Lambda = YukawaModel.Lambda(nion(i)+1,1:YKorder);
        for j = 1:YKorder
            Ffactor(i,:) = Ffactor(i,:) + (Zs-nion(i))*A(j)./(1+(qbar./alpha./Lambda(j)).^2);
        end
    end

    if plotfig==true
        figure
           for i = 1:numel(nion)
            semilogx(qbar, Ffactor(i,:), '-', 'linewidth', 2, 'DisplayName', ['nion=',num2str(nion(i))])
            hold on
           end
            ax = gca;        
            ax.XLim = [1e-4,50];
            ax.YLim = [0, Zs];
            grid on
            grid minor
            set(gca, 'FontSize', 25);
            legend('Location', 'northeast', 'interpreter', 'latex')
            title(['AFF, Z=',num2str(Zs), ', ',num2str(YKorder),' Yukawa exponentials'], 'Interpreter', 'latex');
            xlabel("$\bar{q}$", 'interpreter', 'latex')
            ylabel("$Z\times F_F(\bar{q})$", 'interpreter', 'latex')
    end
    % return AFF and q.
    YukawaModel.qbar = qbar;
    YukawaModel.AFF = Ffactor;
end
