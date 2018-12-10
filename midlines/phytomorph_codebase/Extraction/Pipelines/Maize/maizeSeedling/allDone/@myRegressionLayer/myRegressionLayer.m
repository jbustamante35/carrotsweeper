classdef myRegressionLayer < nnet.layer.RegressionLayer
        
    properties
        % (Optional) Layer properties.

        % Layer properties go here.
    end
 
    methods
        function layer = myRegressionLayer()           
            % (Optional) Create a myRegressionLayer.

            % Layer constructor function goes here.
        end

        function loss = forwardLoss(layer, Y, T)
            % Return the loss between the predictions Y and the 
            % training targets T.
            %
            % Inputs:
            %         layer - Output layer
            %         Y     – Predictions made by network
            %         T     – Training targets
            %
            % Output:
            %         loss  - Loss between Y and T

            % Layer forward loss function goes here.
            
            % Calculate sum of squares.
            T = permute(T,[4 3 1 2]);
            T = reshape(T,size(Y));
            sumSquares = sum((Y-T).^2);
            sumSquares = sum((Y-T));
            % Take mean over mini-batch.
            N = size(Y,1);
            loss = sum(sumSquares(:))/N;
            if loss > 10^4
                here= 1;
            end
             if any(isnan(loss))
                here = 1;
            end
        end
        
        function dLdY = backwardLoss(layer, Y, T)
            % Backward propagate the derivative of the loss function.
            %
            % Inputs:
            %         layer - Output layer
            %         Y     – Predictions made by network
            %         T     – Training targets
            %
            % Output:
            %         dLdY  - Derivative of the loss with respect to the predictions Y        

            % Layer backward loss function goes here.
            T = permute(T,[4 3 1 2]);
            T = reshape(T,size(Y));
            N = size(Y,1);
            dLdY = 2*(Y-T);
            
            if any(isnan(dLdY))
                here = 1;
            end
            dLdY = Y;
        end
    end
end