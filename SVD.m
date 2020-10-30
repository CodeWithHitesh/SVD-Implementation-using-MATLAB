%close all;
clc;
clear;
 
% Reading image as input
I = imread('img/phantom2.tif');

 
f=1;
%neighborhood window size = 2f+1 , ie 5x5
t=3;
% search window size = 2t+1 , ie 11x11
 

% Noisy image.
noisex = imread('img/phantom2_0.2.tif');
noisy = double(noisex);
[m,n] = size(noisex);


% Assign a clear output image
cleared = zeros(m,n);
 

h_1=37;


% Replicate boundaries of noisy image
noisy2 = padarray(noisy,[f,f],'symmetric');



% Now we'll calculate ouput for each pixel
for i=1:m
    for j=1:n
        
        % to compensate for shift due to padarray function
        im = i+f;   
        jn= j+f;
        
        
        % neighborhood of concerned pixel
        W1 = noisy2(im-f:im+f , jn-f:jn+f);
        
        
        
        % Performing SVD over W1.
        [temp1,temp2]=size(W1);
        [EVTU  EVSU]=eig(W1*(W1'));
        evalues_U=round(diag(EVSU));
        
        
        [EVTV  EVSV]=eig((W1')*W1);
        evalues_V=round(diag(EVSV));
        
        
        evalues_U=evalues_U';
        evalues_V=evalues_V';
        
        h=1;
        s1=size(evalues_U,2);
        s2=size(evalues_V,2);
        common=zeros(1,(2*f)+1);
        for i1=1:s1
            for j1=1:s2
                if evalues_U(1,i1)==evalues_V(1,j1)
                    common(1,h)=evalues_V(1,j1);
                    h=h+1;
                    break;
                end
            end
        end
        
        common=sort(common,'descend');
        h=1;
        svdh=zeros(temp1,temp2);
        for i1=1:temp1
            for j1=1:temp2
                if(i1==j1)
                   svdh(i1,j1)=sqrt(common(1,h));
                   h=h+1;
                end
            end
        end
        W1s=diag(svdh);

        
        
        
        
        
        
        % BOundaries of search window for that pixel
        rmin = max(im-t, f+1);
        rmax = min(im+t, m+f);
        smin = max(jn-t, f+1);
        smax = min(jn+t, n+f);
        
        
        NL=0;    
        Z =0;    
        % Run loop through all the pixels in search window
        
        for r=rmin:rmax
            for s=smin:smax
                
                % neighborhood of pixel 'j' being compared with
                % neighbourhood of concerned pixel.
                W2 = noisy2(r-f:r+f, s-f:s+f);
                
                
                
                % Performing svd over W2. 
                [temp1,temp2]=size(W2);
                [EVTU  EVSU]=eig(W2*(W2'));
                evalues_U=round(diag(EVSU));
                
                
                [EVTV  EVSV]=eig((W2')*W2);
                evalues_V=round(diag(EVSV));
                
                
                evalues_U=evalues_U';
                evalues_V=evalues_V';
                
                h=1;
                s1=size(evalues_U,2);
                s2=size(evalues_V,2);
                common=zeros(1,(2*f)+1);
                for i1=1:s1
                    for j1=1:s2
                        if evalues_U(1,i1)==evalues_V(1,j1)
                            common(1,h)=evalues_V(1,j1);
                            h=h+1;
                        end
                    end
                end
                
                common=sort(common,'descend');
                h=1;
                svdh=zeros(temp1,temp2);
                for i1=1:temp1
                    for j1=1:temp2
                        if(i1==j1)
                           svdh(i1,j1)=sqrt(common(1,h));
                           h=h+1;
                        end
                    end
                end
        
                W2s=diag(svdh);

                
                
                
                
                
                % L2 norm.
                d2= sum((W1s-W2s).^2);
                wij=exp(-d2/(h_1*h_1));
                
                % update Z and NL
                Z = Z + wij;
                NL = NL + (wij*noisy2(r,s));
            end
        end
        % normalization of NL
        cleared(i,j) = NL/Z;
    end
end
% convert cleared to uint8
cleared = uint8(cleared);
 

% show results
figure(1);
imshow(noisex),title('noisy Image');

figure(2);
imshow(uint8(I)),title('Orignal Image');

figure(3);
imshow(uint8(cleared)),title('output of NLM-SVD');


psn1=PS(noisex,I);
psn2=PS(cleared,I);

%Program for Peak Signal to Noise Ratio Calculation
function PSR = PS(I, D)

I = double(I);
D = double(D);

[M N] = size(I);
error = I - D;
MSE = sum(sum(error .* error)) / (M * N);
MSE
if(MSE > 0)
    PSR = 10*log(255*255/MSE) / log(10);
else
    PSR = 99;
end
end
