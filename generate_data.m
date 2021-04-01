%%clear
clear;
tic

%% initialize values
sigma = 2.5;                    %filtering strength
vertice = 32;                   %vertice length of cube
level = 0.517;                  %cut-off level that determines porosity
boundaries = [10,35];           %min and max porosity
n = 100000;                     %amount of porous images
i = 0;                          %position of stored data
iteration = 0;                  %counting iterations per successful geometry

%% generation
while 1
    image = rand(vertice,vertice,vertice); %generate random noise in a 3D matrix
    filtered = imgaussfilt3(image,sigma);  %apply a gaussian averaging filter
    BW = imbinarize(filtered,level);       %cut off and replace with binary values
    
    %Determine if material is porous in z-direction
    locations = find(BW(:,:,1));
    BW2 = imfill(~BW,locations);
    BW3 =  BW - ~BW2;  %The change after z-direction filling algorithm
    
    %Determine if material is porous in x-direction
    locations2 = find(BW(:,1,:));
    locations2 = floor((locations2-1)/vertice)*vertice^2 + 1 + mod(locations2-1,vertice);
    BW22 = imfill(~BW,locations2);
    BW32 =  BW - ~BW22; %The change after x-direction filling algorithm
    
    %Determine if material is porous in y-direction
    locations3 = find(BW(1,:,:));
    locations3 = floor((locations3-1)/vertice)*vertice^2 + 1 + vertice*mod(locations3-1,vertice);
    BW23 = imfill(~BW,locations3);
    BW33 =  BW - ~BW23; %The change after y-direction filling algorithm
    
    iteration = iteration + 1;
    
    %determine percentage of ones / porosity
    percentageOfOnes = sum(BW(:) == 1) / numel(BW) * 100;
    
    %boundary conditions in the if statement: porous in all 3 directions,
    %and the proper porosity value.
    if (sum(BW3(:,:,vertice),'all') ~= 0 && sum(BW32(:,vertice,:),'all') ~= 0 && sum(BW33(vertice,:,:),'all') ~= 0 && percentageOfOnes > boundaries(1) && percentageOfOnes < boundaries(2))
        i = i + 1;
        dataset(:,:,:,i) = BW;              %store data
        level = unifrnd(0.506,0.517);       %refresh the cut-off boundary for a different porisity
        iterations(1,i) = iteration;
        percentages(i,1) = percentageOfOnes;
        iteration = 0;
        if i == n
            break                           %quit the program when the dataset is full
        end
        
    end
    
end
toc

save('dataset.mat', 'dataset','percentages','-v7.3');

%% plotting
% Uncomment these lines for plots.

% figure(1)
% colormap gray
% slice(double(BW),1:vertice,1:vertice,1:vertice)
%
% figure(2)
% colormap gray
% slice(double(filtered),1:vertice,1:vertice,1:vertice)
%
% figure(3)
% colormap gray
% slice(double(image),1:vertice,1:vertice,1:vertice)
% 
% figure(4)
% histogram(percentages,30)

