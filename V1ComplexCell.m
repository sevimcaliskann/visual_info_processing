classdef V1ComplexCell <handle
    properties(Access = 'public')
        orientation, 
        firstDeriv,
        secondDeriv,
        responseFrames
    end
    
    
    methods(Access='public')
        function this = V1ComplexCell(matrixSize, scale, orientation, frequency)
            this.orientation = orientation;
            pointX = ceil(matrixSize(1)/2); 
            pointY = ceil(matrixSize(2)/2);
          
            maxX = matrixSize - pointX;
            maxY = matrixSize - pointY;
          
            if(mod(matrixSize, 2) ==0) %if matrix size is even
                maxX = maxX - 1; %kernel size should be odd
                maxY = maxY - 1;
            end

            [x, y] = meshgrid(-pointX :1: maxX, -pointY :1: maxY);
%             x = double(x); y = double(y);
            this.firstDeriv = this.differentiation(scale, orientation, frequency,1,x,y);
            this.secondDeriv = this.differentiation(scale, orientation, frequency,2,x,y);
            
        end
        
        function ShowResponse(this)			
            figure('Name', 'Response');
            if size(this.responseFrames,3)>0
                imshow(this.responseFrames(:,:,1), []);
            end
        end
   
        
        function ShowFilter(this, featureExtractionFunction, displayFunction)
            
			if ~exist('featureExtractionFunction', 'var') || isempty(featureExtractionFunction)
				featureExtractionFunction = @(x) GaborKernel.GetRealParts(x);
			end
			
			if ~exist('displayFunction', 'var') || isempty(displayFunction)
				displayFunction = @(im) surf(im);
			end
			
            figure('Name', 'Gabor Kernels');
                
            filterValues = this.firstDeriv; 
            filterValues = featureExtractionFunction(filterValues);
            filterValues = imresize(filterValues, [30 30]);
            subplot(1, 2, 1)                    
            displayFunction(filterValues);
            az = 135;
            el = 60;
            view(az, el);
            
            
            filterValues = this.secondDeriv; 
            filterValues = featureExtractionFunction(filterValues);
            filterValues = imresize(filterValues, [30 30]);
            subplot(1, 2, 2)                    
            displayFunction(filterValues);
            az = 135;
            el = 60;
            view(az, el);
            
        end
        
        
        function kernel = differentiation(this, scale, orientation, frequency, degree, xCoord, yCoord)        
           syms x y
           e = exp(1); 
            
           r = 1; %constant (values ~=1 -> eliptical filter)

           R1 = x.*cos(orientation) + y.*sin(orientation);
           R2 =-x.*sin(orientation) + y.*cos(orientation);

           expFactor = -1/2 * ( (R1/scale).^2 + (R2/(r*scale)).^2  );

           gauss = 1 / ( sqrt(r*pi)*scale) ;
           gauss =  gauss .* e.^expFactor;

           gaborReal = gauss .* cos(frequency*R1);
           gaborImag = gauss .* sin(frequency*R1);
           for i = 1:degree        
               gaborReal = diff(gaborReal, x, y);
               gaborImag = diff(gaborImag, x, y);
           end
           gaborReal = matlabFunction(gaborReal);
           gaborImag = matlabFunction(gaborImag);
           
           
           
           kernel = gaborReal(xCoord, yCoord) + gaborImag(xCoord, yCoord)*1i;
           kernel = this.RemoveDC(kernel);
           kernel = kernel'; %we want kernel(y,x); y=row; x=column;
        end
        
        
        function rez = RemoveDC(this, matrix)
            
            valueDC = mean( matrix(:) );
            rez = matrix - valueDC;
        end
        
        
        function out = Integrate(this, frames, duration, featureExtractionFunction)
            if ~exist('featureExtractionFunction', 'var') || isempty(featureExtractionFunction)
				featureExtractionFunction = @(x)GaborKernel.GetAmplitudes(x);
            end
            
            
            this.responseFrames = [];
            for i = 1:duration:size(frames,3)
                responseFrame = this.IntegrateAlongDuration(frames, i, min(i+duration-1, size(frames,3)));
                responseFrame = featureExtractionFunction(responseFrame);
                this.responseFrames = cat(3, this.responseFrames, responseFrame);
            end
            out = this.responseFrames;
        end
        
        function out = IntegrateAlongDuration(this, frames, start, ending)
            tauv = ending - start +1;
            rFrames = [];
            for i = start:ending
                filt = temporalFilter(this, i, tauv);
                rFrame = conv2(frames(:,:,i), filt, 'same');
                rFrames = cat(3, rFrames, rFrame);
            end
            
            out = sqrt(sum(rFrames.^2, 3));
        end
        
        function out = gamma(this, t, tauv, n)
            out = ((t^n)/(tauv^(n+1)*factorial(n)))*exp(-t/tauv);
        end
        
        function out = temporalHFunction(this, typeStr, t, tauv) % t = time(frame number), tav = integration duration
            if strcmp(typeStr, 'fast')==1
                out = gamma(this, t, tauv, 3) - gamma(this, t, tauv, 5);
            elseif strcmp(typeStr, 'slow') ==1 
                out = gamma(this, t, tauv, 5) - gamma(this, t, tauv, 7);
            else
                out = 0;
            end
        end
        
        function out = temporalFilter(this, t, tauv)

            hfast = this.temporalHFunction('fast', t, tauv);
            hslow = this.temporalHFunction('slow', t, tauv);
            
            out = hfast*this.firstDeriv - hslow*this.secondDeriv;
            
        end
        
        
    end
    
    
end