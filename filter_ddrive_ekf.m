% Localization algorithm for a Differential Drive Mobile Robot based
% on the Extended Kalman Filter (EKF)

function filter = filter_ddrive_ekf()
    filter = filter_localization2d(@filterStep); % reuse drawing function from generic localization2d block     
    filter.depends = {'sensors/odometer', 'sensors/landmark_detector', 'environment/landmarks'};
    
    filter.default_initialPose = [0 0 0]';
    filter.default_initialPoseCov = zeros(3, 3);
    
    filter.default_odometryError = 5 * pi / 180;    % sigma of assumed odometer uncertainty in rad/s
    filter.default_wheelRadiusError = 1e-3;         % sigma of assumed wheel diameter uncertainty in m
    filter.default_wheelDistanceError = 5e-3;       % sigma of assumed wheel distance uncertainty in m
    filter.default_bearingError = 5 * pi / 180;     % sigma of assumed bearing error in rad
    filter.default_rangeError = 2 / 100;            % sigma of assumed range error in percent of the aquired distance value
    
end

function [state, out, debugOut] = filterStep(block, t, state, odometer, sensor, landmarks)
    debugOut = [];

    if isempty(state)            
        % Initialization
        state.pose = block.initialPose(:);
        state.cov = block.initialPoseCov;
        state.lastInput = [0; 0]; % initial speed is zero
        state.t = 0;
    end

    % use shorter symbols 
    x = state.pose;
    P = state.cov;
    u = state.lastInput;
    tNow = state.t;

    iPredict = 1;
    iUpdate = 1;

    while tNow < t
        % determine, which measurement to proceed next
        if iUpdate <= length(sensor)
            tNextUpdate = sensor(iUpdate).t;
        else tNextUpdate = t;
        end            
        while tNow < tNextUpdate
            if iPredict <= length(odometer)
                tNext = odometer(iPredict).t;
                if tNext <= tNextUpdate
                    [x, P] = doPrediction(x, P, u, tNext - tNow);
                    tNow = tNext;           
                    u = odometer(iPredict).data(:);
                    iPredict = iPredict + 1;
                else break;
                end
            else break;    
            end
        end

        if tNow < tNextUpdate
            [x, P] = doPrediction(x, P, u, tNextUpdate - tNow);
            tNow = tNextUpdate;
        end

        if iUpdate <= length(sensor)                
            [x, P] = doUpdate(x, P, sensor(iUpdate).data, landmarks.data);
            iUpdate = iUpdate + 1;
        end
    end

    % put short-named intermediate variables back into the filter state
    state.pose = x;
    state.cov = P;
    state.lastInput = u;
    state.t = tNow;

    % the output of the localization filter is the estimated state (=pose) vector 
    out.pose = x;  
    out.cov = P;

    function [x, P] = doPrediction(x, P, u, T)
        % Implementation of the prediction step
      
        % get the model parameters 
        if numel(block.wheelRadius) == 1
            R_R = block.wheelRadius;
            R_L = block.wheelRadius;
        else
            R_R = block.wheelRadius(1);
            R_L = block.wheelRadius(2);
        end
        a = block.wheelDistance / 2;
        
        dtheta_R = u(1);
        dtheta_L = u(2);
        
        % TODO: implement the prediction step
        v = (R_R * dtheta_R + R_L * dtheta_L) / 2;
        w = (R_R * dtheta_R - R_L * dtheta_L) /(2 * a);
        
%         A = [0 0 -sin(x(3))*v
%              0 0  cos(x(3))*v
%              0 0            0];
%         
%         B = [R_R*cos(x(3))/2, R_L*cos(x(3))/2
%              R_R*sin(x(3))/2, R_L*sin(x(3))/2
%              R_R/2/a        ,-R_L/2/a        ];
         
        F = [1 0 -v*T*sin(x(3))
             0 1  v*T*cos(x(3))
             0 0              1];
         
        H = [ (R_R*T*cos(x(3)))/2 - (R_R*T^2*v*sin(x(3)))/(4*a), (R_L*v*sin(x(3))*T^2)/(4*a) + (R_L*cos(x(3))*T)/2
              (R_R*v*cos(x(3))*T^2)/(4*a) + (R_R*sin(x(3))*T)/2, (R_L*T*sin(x(3)))/2 - (R_L*T^2*v*cos(x(3)))/(4*a)
              (R_R*T)/(2*a)                                    ,                                    -(R_L*T)/(2*a)];
          
%         if (w ==0)
%             F = [1 0 -v*T*sin(x(3))
%                  0 1  v*T*cos(x(3))
%                  0 0              1];
%               
%             H = [R_R*cos(x(3))*T/2, R_L*cos(x(3))*T/2
%                  R_R*sin(x(3))*T/2, R_L*sin(x(3))*T/2
%                  R_R*T/2/a          ,-R_L*T/2/a          ];
%         else
%             F = [1 0 (cos(w*T + x(3))-cos(x(3)))/w*v
%                  0 1 (sin(w*T + x(3))-sin(x(3)))/w*v
%                  0 0                               1];
%                        
%             H = [((T*v*cos(x(3) + T*w))/w - (v*(sin(x(3) + T*w) - sin(x(3))))/w^2) * R_R/2/a + (sin(x(3) + T*w) - sin(x(3)))*R_R/2/w, -((T*v*cos(x(3) + T*w))/w - (v*(sin(x(3) + T*w) - sin(x(3))))/w^2) * R_L/2/a + (sin(x(3) + T*w) - sin(x(3)))*R_L/2/w
%                  ((v*(cos(x(3) + T*w) - cos(x(3))))/w^2 + (T*v*sin(x(3) + T*w))/w) * R_R/2/a - (cos(x(3) + T*w) - cos(x(3)))*R_R/2/w, -((v*(cos(x(3) + T*w) - cos(x(3))))/w^2 + (T*v*sin(x(3) + T*w))/w) * R_L/2/a - (cos(x(3) + T*w) - cos(x(3)))*R_L/2/w
%                  R_R*T/2/a                                                                                                            , -R_L*T/2/a                                                                                                           ];           
%         end
        if block.useNumericPrediction
            options = odeset('RelTol',1e-13,'AbsTol',[1e-4 1e-4 1e-5]);
            [~,X] = ode45(@(t,x)kinematicmodel(t,x,v,w),[0 T],x,options);
            x = X(end,:).';
        end
        
        if block.useExactDiscretization
            if (w == 0)
                x(1) = x(1) + v*T*cos(x(3));
                x(2) = x(2) + v*T*sin(x(3));
            else
                x(1) = x(1) + ( sin(w*T + x(3)) - sin(x(3)))/w*v;
                x(2) = x(2) + (-cos(w*T + x(3)) + cos(x(3)))/w*v;
            end      
            x(3) = x(3) + w*T;
        end
        if x(3) >  pi
            x(3) = x(3) - 2*pi;
        end
        if x(3) < -pi
            x(3) = x(3) + 2*pi;
        end
        N =blkdiag((block.odometryError)^2, (block.odometryError)^2);
        P = F*P*F' + H*N*H'; 
        
        function dx = kinematicmodel(t,x,v,w)
            dx = zeros(3,1);
            dx(1) =  v*cos(x(3));
            dx(2) =  v*sin(x(3));
            dx(3) =  w;
        end
    end

    function [x, P] = doUpdate(x, P, meas, landmarks)
        % Implementation of the update step
        visIdx = meas.lmIds; % assume we can associate each measurement with a known landmark
        if isempty(visIdx); return; end
        
        beta = meas.bearing;
        d = meas.range;                        
        m = landmarks(visIdx, :);

        % TODO: implement the update step
        x1 = sym('x1');
        x2 = sym('x2');
        x3 = sym('x3');
        
        d_tilt = sqrt((m(:,1) - x1).^2 + (m(:,2) - x2).^2);
        beta_tilt = atan2(m(:,2) - x2, m(:,1) - x1) - x3;
        C1 = jacobian(d_tilt   , [x1, x2, x3]);
        C2 = jacobian(beta_tilt, [x1, x2, x3]);
        C1 = eval(subs(C1, 'x1,x2,x3',{x(1),x(2),x(3)}));
        C2 = eval(subs(C2, 'x1,x2,x3',{x(1),x(2),x(3)}));
        
        d_tilt    = eval(subs(d_tilt   , 'x1,x2,x3',{x(1),x(2),x(3)}));        
        beta_tilt = eval(subs(beta_tilt, 'x1,x2,x3',{x(1),x(2),x(3)}));
        temp1 = beta_tilt >  pi;
        temp2 = beta_tilt < -pi;
        beta_tilt = beta_tilt - temp1*2*pi + temp2*2*pi;
        W1 =diag(block.rangeError  *block.rangeError  *ones(1,length(visIdx)));
        W2 =diag(block.bearingError*block.bearingError*ones(1,length(visIdx)));

        Error1 = d - d_tilt;
        Error2 = beta - beta_tilt;

        if (block.useRange && block.useBearing)
            C=[C1;C2];
            W = blkdiag(W1,W2);
            Error = [Error1; Error2];
        elseif (block.useRange)
            C = C1;
            W = W1;
            Error = Error1;
        elseif (block.useBearing)
            C = C2;
            W = W2;
            Error = Error2;
            save C.mat C
        end
        
        K = P*C'/(C*P*C' + W); 
        x = x + K*Error;
        P = (eye(3) - K*C) * P;
    end
end

