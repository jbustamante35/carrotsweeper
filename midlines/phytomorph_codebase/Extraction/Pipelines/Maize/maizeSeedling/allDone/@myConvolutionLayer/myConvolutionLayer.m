classdef myConvolutionLayer < Convolution2D
    
    
    
    
    
    methods
        
        function this = myConvolutionLayer(varargin)
            this@Convolution2D(varargin{:})
        end
        
        function this = set.Weights(this, weights)
            this.LearnableParameters(this.WeightsIndex) = weights;
            if ~this.CacheHandle.isEmpty
                this.CacheHandle.clearCache;
            end
        end
    end
end