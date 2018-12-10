classdef locationRegressionLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties

        % Layer properties go here
    end

    properties (Learnable)
        % (Optional) Layer learnable parameters
        kernelStack;
        % Layer learnable parameters go here
    end
    
    methods
        function layer = locationRegressionLayer()
            
        end
        
        function Z = predict(layer, X)
            Z = layer.forward(X);
        end

        function [Z, memory] = forward(layer, X)
            for batch = 1:size(X,4)
                for channel = 1:size(3)
                    Z(:,:,channel,batch) = conv2(X,layer.kernelStack));
                end
            end
        end

        function [dLdX, dLdW1, ..., dLdWn] = backward(layer, X, Z, dLdZ, memory)
                
            % Backward propagate the derivative of the loss function through 
            % the layer
            %
            % Inputs:
            %         layer             - Layer to backward propagate through
            %         X                 - Input data
            %         Z                 - Output of layer forward function            
            %         dLdZ              - Gradient propagated from the deeper layer
            %         memory            - Memory value which can be used in
            %                             backward propagation
            % Output:
            %         dLdX              - Derivative of the loss with respect to the
            %                             input data
            %         dLdW1, ..., dLdWn - Derivatives of the loss with respect to each
            %                             learnable parameter
            
            % Layer backward function goes here
            
        end
    end
end