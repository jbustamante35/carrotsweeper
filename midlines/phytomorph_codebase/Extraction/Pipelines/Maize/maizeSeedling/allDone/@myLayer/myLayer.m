classdef myLayer < nnet.layer.Layer

    properties
        % (Optional) Layer properties.
        h;
        % Layer properties go here.
    end

    properties (Learnable)
        % (Optional) Layer learnable parameters.

        % Layer learnable parameters go here.
    end
    
    methods
        function layer = myLayer(h)
            % (Optional) Create a myLayer.
            % This function must have the same name as the layer.
            layer.h = h;
            % Layer constructor function goes here.
        end
        
        function Z = predict(layer, X)
            % Forward input data through the layer at prediction time and
            % output the result.
            %
            % Inputs:
            %         layer    -    Layer to forward propagate through
            %         X        -    Input data
            % Output:
            %         Z        -    Output of layer forward function
            
            % Layer forward function for prediction goes here.
             here = 1;
        end

        function [Z, memory] = forward(layer, X)
            % (Optional) Forward input data through the layer at training
            % time and output the result and a memory value.
            %
            % Inputs:
            %         layer  - Layer to forward propagate through
            %         X      - Input data
            % Outputs:
            %         Z      - Output of layer forward function
            %         memory - Memory value for backward propagation

            memory = [];
            
            % Layer forward function for training goes here.
            [p2,p1] = ndgrid(1:size(X,1),1:size(X,2));
            
            for tr = 1:size(X,4)
                for k = 1:size(X,3)
                    tmp = (X(:,:,k,tr));
                    %tmp = .5*(1+tmp);
                    h(k,tr) = sum(tmp(:));
                    %{
                    h(k,tr);
                    if h(k,tr) == 0
                        
                        h(k,tr) = rand(1);
                        tmp = rand(size(tmp));
                        h(k,tr) = sum(tmp(:));
                        zero = 1;
                    end
                    %}
                    if ~isempty(layer.h)
                        if tr == 1
                            figure(layer.h(1));
                            imshow(tmp,[]);
                            drawnow
                        end
                    end
                    
                    g1(k,tr) = p1(:)'*tmp(:);
                    g2(k,tr) = p2(:)'*tmp(:);
                    Z(tr,k,1) = g1(k,tr).*h(k,tr).^-1;
                    %Z(tr,k,2) = g2(k,tr).*h(k,tr).^-1;
                    dZdX(1,:,:,k,tr) = (p1*h(k,tr) - g1(k,tr))*h(k,tr)^-2;
                    %dZdX(2,:,:,k,tr) = (p2*h(k,tr) - g2(k,tr))*h(k,tr)^-2;
                    %{
                    tmp = tmp / sum(tmp(:));
                    tmp = imfilter(tmp,fspecial('gaussian',[11 11],5),'replicate');
                    
                    
                    [dtmpd1,dtmpd2] = gradient(tmp);
                    
                    idtmpd1 = dtmpd1.^-1;
                    idtmpd2 = dtmpd2.^-1;
                    TV = .2;
                    idtmpd1(idtmpd1 > TV) = TV;
                    idtmpd2(idtmpd2 > TV) = TV;
                    idtmpd1(idtmpd1 < -TV) = -TV;
                    idtmpd2(idtmpd2 < -TV) = -TV;
                    idtmpd1 = imfilter(idtmpd1,fspecial('gaussian',[11 11],5),'replicate');
                    idtmpd2 = imfilter(idtmpd2,fspecial('gaussian',[11 11],5),'replicate');
                    
                    dZdX(1,:,:,k,tr) = p1.*idtmpd1;
                    dZdX(2,:,:,k,tr) = p2.*idtmpd2;
                    %}
                end
            end
            dZdX(isnan(dZdX)) = 0;
            dZdX(isinf(dZdX)) = 0;
            memory = dZdX;
            if any(isnan(dZdX(:)))
                here = 1;
            end
            if any(Z(:) > 20^4)
                here = 1;
            end
        end

        function [dLdX] = backward(layer, X, Z, dLdZ, memory)
            % Backward propagate the derivative of the loss function through 
            % the layer.
            %
            % Inputs:
            %         layer             - Layer to backward propagate through
            %         X                 - Input data
            %         Z                 - Output of layer forward function            
            %         dLdZ              - Gradient propagated from the deeper layer
            %         memory            - Memory value from forward function
            % Outputs:
            %         dLdX              - Derivative of the loss with respect to the
            %                             input data
            %         dLdW1, ..., dLdWn - Derivatives of the loss with respect to each
            %                             learnable parameter
            try
                % Layer backward function goes here.
                for tr = 1:size(X,4)
                    for k = 1:size(X,3)
                        dLdX(:,:,k,tr) = memory(1,:,:,k,tr).*dLdZ(tr,k,1);% + memory(2,:,:,k,tr).*dLdZ(tr,k,2);
                        %dLdX(:,k,tr) = memory(:,:,k,tr)'*squeeze(dLdZ(tr,k,:));
                    end




                end
                %dLdX = permute(dLdX,[3 2 1]);
                %sz = size(X);
                %dLdX = reshape(dLdX,sz);
                %dLdX = dLdX * 100000;

                if ~isempty(layer.h)
                    figure(layer.h(2));
                    imshow(dLdX(:,:,1),[]);
                    drawnow
                end



                if any(isnan(dLdX(:)))
                    here = 1;
                end
            catch ME
                ME;
            end
        end
    end
end