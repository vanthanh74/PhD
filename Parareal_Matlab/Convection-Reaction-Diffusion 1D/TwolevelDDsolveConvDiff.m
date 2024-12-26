function [U,L2NormError,LInfNormError] = TwolevelDDsolveConvDiff(f,A,U_fine,N,dX,dT,R1,A1,R_sub,A_sub,R,B,pick,nx_fine,uF_seq_mat,type)



% FC - relaxation
% residual
r_k = f - A*U_fine;

% first level ~ fine propagator ~ additive Schwarz in parallel
%     U_fine_temp = P\r_k;
for i = 1:N+1
    %         i = i + 1
    if i == 1
        x_temp_1 = R1*r_k;
        x_temp_2 = A1\x_temp_1;
        x{i} = R1'*x_temp_2;
        U_fine_temp = x{i};
    else
        x_temp_1 = R_sub{i}*r_k;
        x_temp_2 = A_sub{i}\x_temp_1;
        x{i} = R_sub{i}'*x_temp_2;
        U_fine_temp = U_fine_temp + x{i};
    end
end

% second level ~ coarse grid correction sequentially
U_fine_temp_1 = R*U_fine_temp;
U_fine_temp_2 = B\U_fine_temp_1;
U_fine_temp_3 = R'*U_fine_temp_2;
U_fine_temp_4 = U_fine_temp_3 + U_fine_temp - R'*R*U_fine_temp;
U_fine = U_fine + U_fine_temp_4;

if strcmp(type,'FC') == 1
    U = U_fine;
elseif strcmp(type,'FCF') == 1
    % %     % residual
    r_k = f - A*U_fine;
    % %     first level ~ fine propagator ~ additive Schwarz in parallel
    % %         U_fine_temp = P\r_k;
    for i = 1:N+1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp1 = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp1 = U_fine_temp1 + x{i};
        end
    end
    U_fine = U_fine + U_fine_temp1;
elseif strcmp(type,'FCFF') == 1
    % % %         residual
    r_k = f - A*U_fine;
    % %     first level ~ fine propagator ~ additive Schwarz in parallel
    % %         U_fine_temp = P\r_k;
    for i = 1:N+1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp2 = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp2 = U_fine_temp2 + x{i};
        end
    end
    U_fine = U_fine + U_fine_temp2;
    
    % % %         residual
    r_k = f - A*U_fine;
    % %     first level ~ fine propagator ~ additive Schwarz in parallel
    % %         U_fine_temp = P\r_k;
    for i = 1:N+1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp2 = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp2 = U_fine_temp2 + x{i};
        end
    end
    U_fine = U_fine + U_fine_temp2;
elseif strcmp(type,'FCFC') == 1
    
    % %      %FCFC
    % residual
    r_k = f - A*U_fine;
    
    %   % first level ~ fine propagator ~ additive Schwarz in parallel
    % %     U_fine_temp = P\r_k;
    for i = 1:N+1
        % %         i = i + 1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp = U_fine_temp + x{i};
        end
    end
    
    % %   % second level ~ coarse grid correction sequentially
    U_fine_temp_1 = R*U_fine_temp;
    U_fine_temp_2 = B\U_fine_temp_1;
    U_fine_temp_3 = R'*U_fine_temp_2;
    U_fine_temp_4 = U_fine_temp_3 + U_fine_temp - R'*R*U_fine_temp;
    U_fine = U_fine + U_fine_temp_4;
    
elseif strcmp(type,'FCFCF') == 1
    % %      %FCFC
    % residual
    r_k = f - A*U_fine;
    
    %   % first level ~ fine propagator ~ additive Schwarz in parallel
    % %     U_fine_temp = P\r_k;
    for i = 1:N+1
        % %         i = i + 1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp = U_fine_temp + x{i};
        end
    end
    
    % %   % second level ~ coarse grid correction sequentially
    U_fine_temp_1 = R*U_fine_temp;
    U_fine_temp_2 = B\U_fine_temp_1;
    U_fine_temp_3 = R'*U_fine_temp_2;
    U_fine_temp_4 = U_fine_temp_3 + U_fine_temp - R'*R*U_fine_temp;
    U_fine = U_fine + U_fine_temp_4;
    
    % %         % residual
    r_k = f - A*U_fine;
    % %     % first level ~ fine propagator ~ additive Schwarz in parallel
    % %     %     U_fine_temp = P\r_k;
    for i = 1:N+1
        if i == 1
            x_temp_1 = R1*r_k;
            x_temp_2 = A1\x_temp_1;
            x{i} = R1'*x_temp_2;
            U_fine_temp3 = x{i};
        else
            x_temp_1 = R_sub{i}*r_k;
            x_temp_2 = A_sub{i}\x_temp_1;
            x{i} = R_sub{i}'*x_temp_2;
            U_fine_temp3 = U_fine_temp3 + x{i};
        end
    end
    U_fine = U_fine + U_fine_temp3;
end


U_fine_mat = [];
for i = 1:length(U_fine)/(nx_fine-1)
    U_fine_mat(:,i) = [U_fine((i-1)*(nx_fine-1)+1:(i)*(nx_fine-1))];
end
U_fine_mat_reduce =  U_fine_mat(:,pick);
U = U_fine;

L2NormError = [sqrt(sum(dX*sum(dT*(U_fine_mat_reduce' - uF_seq_mat').^2)));];%sqrt(sum(sum(dT*(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)).^2)));%
LInfNormError = [norm(U_fine_mat_reduce' - uF_seq_mat',inf)];%max(abs(uF_seq(end-nx_coarse+2:end,:)-U(end-nx_coarse+2:end,:)));%














