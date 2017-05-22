classdef V1SimpleCell <handle
    properties(Access = 'public')
        
        
        matrixSize,
        
		scale=1,
		orientation,
        frequency,
        
        centerPoint, 
        responseImage
    end
    
    properties(Access = 'public')
        Filter,
        FilterParams,
		FilterValues
    end
    
    
    methods(Access='public')
        function this = V1SimpleCell(matrixSize, orientation, frequency)
            this.matrixSize = matrixSize;
            this.orientation = orientation;
            this.frequency = frequency;
            this.centerPoint = [ceil(matrixSize/2) ceil(matrixSize/2)];
            this.Filter = GaborKernel(this.matrixSize, this.scale, this.orientation, this.frequency, this.centerPoint);
            filterParams = struct('CenterPoint', this.centerPoint, ...
                                      'Scale', this.scale, 'Orientation', this.orientation, ...
									  'Frequency', this.frequency);
									  
            this.FilterParams = filterParams;						
            this.FilterValues = (this.Filter.KernelValues);
        end
        
        function ShowResponse(this)			
            figure('Name', 'Response');    
            imshow(this.responseImage, []);
            
        end
   
        
        function ShowFilter(this, featureExtractionFunction, displayFunction)
            
			if ~exist('featureExtractionFunction', 'var') || isempty(featureExtractionFunction)
				featureExtractionFunction = @(x) GaborKernel.GetRealParts(x);
			end
			
			if ~exist('displayFunction', 'var') || isempty(displayFunction)
				displayFunction = @(im) surf(im);
			end
			
            figure('Name', 'Gabor Kernel');
            filterParams = this.FilterParams;
                
            filterValues = this.FilterValues; 
            filterValues = featureExtractionFunction(filterValues);
            filterValues = imresize(filterValues, [30 30]);
                                
            displayFunction(filterValues);
            az = 135;
            el = 60;
            view(az, el);
            title(strcat( num2str(filterParams.Scale), {'; '}, num2str(filterParams.Frequency), {'; '}, num2str(filterParams.Orientation) ));
            
        end
        
        function response = Convolve(this, image, featureExtractionFunction)
		    if ~exist('featureExtractionFunction', 'var') || isempty(featureExtractionFunction)
				featureExtractionFunction = @(x)GaborKernel.GetAmplitudes(x);
            end
				
            filterValues = this.FilterValues; 
			 						
			response = conv2(image, filterValues, 'valid');
			response = featureExtractionFunction(response);
            this.responseImage = response;			
		end
        
        function out = gamma(t, tauv, n)
            out = ((t^n)/(tauv^(n+1)*factorial(n)))*exp(-t/tauv);
        end
        
        function out = temporalHFunction(typeStr, t, tauv)
            if strcmp(typeStr, 'fast')==1
                out = gamma(t, tauv, 3) - gamma(t, tauv, 5);
            elseif strcmp(typeStr, 'slow') ==1 
                out = gamma(t, tauv, 5) - gamma(t, tauv, 7);
            else
                out = 0;
            end
        end
        
        function out = temporalResponse(this)
            
            g = gabor(this.frequency,this.orientation);

            out = diff(g);
        end
        
        
    end
    
    
end