%DISCLAIMER: This code contains functions that do not have a
%CUDA-equivalent. In this state, the code cannot be ran on a GPU using for
%example the GPU coder, but it does use the parallel processing toolbox of
%matlab for multicore CPU processing.

clear;
format long g

%% initialize values
load('dataset.mat');                  %Load 32x32x32xN 3D matrices and porosity values from 'generate_data.m'
vertice = length(dataset(:,:,:,1));   %Find the length of the cube's vertices
n = length(dataset(1,1,1,:));         %Find the amount of cubes in the set
depthnumber = vertice^3 - vertice^2;  %depthnumber and verticesq are the same for the entire set, only calculate once
verticesq = vertice^2;
pressure = [vertice+1,0];             %solutions independent of pressure, this gives a gradient of 1 anyhow.

perm = zeros(8*n,1);                  %predefining the permeability vector.

tic

%% Boundary conditions and finding k
parfor v=1:n                          %Parfor enables hyperthreading (use the parallel processing toolbox).
    for xyz=1:3                       %3 loops for calculations in 3 directions
        geometry = dataset(:,:,:,v);
        
        %Here the locations of the x,y, and z planes are found. They are
        %automatically mapped to the 1-1024 plane so the second floor()
        %function maps them back to their actual position.
        if xyz==1
            locations = find(geometry(:,1,:));  
            locations = floor((locations-1)/vertice)*vertice^2 + 1 + mod(locations-1,vertice);
            locations2 = find(geometry(:,vertice,:)); 
            locations2 = floor((locations2-1)/vertice)*vertice^2 + vertice*(vertice-1) + 1 + mod((locations2-1),vertice);
        elseif xyz==2
            locations = find(geometry(1,:,:)); 
            locations = floor((locations-1)/vertice)*vertice^2 + 1 + vertice*mod((locations-1),vertice);
            locations2 = find(geometry(vertice,:,:)); 
            locations2 = floor((locations2-1)/vertice)*vertice^2 + vertice + vertice*mod((locations2-1),vertice);
        elseif xyz==3
            locations = find(geometry(:,:,1));
            locations2 = find(geometry(:,:,vertice))+depthnumber;
        end
        %Fillforward and fillbackward together determine the unreachable or
        %unused voxels of the structure, these are removed to decrease
        %computation time as they do not contribute.
        fillforw = imfill(~geometry,locations);
        fillback = imfill(~geometry,locations2);
        unreachable = ~fillforw | ~fillback;
        geometry=  geometry - unreachable;
        
        
        %position finds the indices of open spaces (1's) and shrinks the
        %size by 65%-90% to decrease computation time.
        position = find(geometry);
        sparsematrix = zeros(length(position),length(position));   %predefining the matrix that contains linear equations of 2-step method.
        keffmatrix = zeros(length(position),length(position));     %predefining the matrix that contains linked permeabilities from voxel to voxel.
        neighbourtype = zeros(length(position),6);                 %for each voxels this matrix will contain it's 3D neighbours.
        
        
        % edge-planes
        yplanetop = [];
        yplanebot = [];
        xplanetop = [];
        xplanebot = [];
        zplanetop = [];
        zplanebot = [];
        
        %This loop runs for every voxel that is open.
        for j=1:length(position)
            %Matlab matrix counting from find() goes as follows;
            %(1,1,1) -> (1,2,1) -> (2,1,1) -> (2,2,1) -> (1,1,2)
            %initialize these values once so they will not be computed multiple
            %times
            mody = mod(position(j),vertice);
            modx = mod(position(j),verticesq);
            %set the diagonal that will be reduced if at edges
            diagonal = 6;
            %y-axis check in positive direction
            if mody ~= 0 %modulo determines if at y-bottom
                if ismember(position(j)+1,position)==1 %is there something under j
                    sparsematrix(j,find(position==position(j)+1)) = -1;
                    neighbourtype(j,1) = j+1;
                end
            else
                diagonal = diagonal -1; %if j is at y-bottom diagonal is decreased
                yplanebot(end+1) = j;
            end
            %y-axis check in negative direction
            if mody ~= 1 %modulo determines if at y-top
                if ismember(position(j)-1,position)==1  %is there something above j
                    sparsematrix(j,find(position==position(j)-1)) = -1;
                    neighbourtype(j,2) = j-1;
                end
            else
                diagonal = diagonal -1; %if j is at y-top diagonal is decreased
                yplanetop(end+1) = j;
            end
            %x-axis check in positive direction
            if modx <= (vertice^2 - vertice) && modx ~= 0 %modulo determines if at x top
                if ismember(position(j)+vertice,position)==1  %is there something to the right of j
                    sparsematrix(j,find(position==position(j)+vertice)) = -1;
                    neighbourtype(j,3) = find(position==position(j)+vertice);
                end
            else
                diagonal = diagonal -1; %if j is at x-edge diagonal is decreased
                xplanetop(end+1) = j;
            end
            %x-axis check in negative direction
            if modx > vertice || modx == 0 %modulo determines if at x bottom
                if ismember(position(j)-vertice,position)==1 %is there something to the left of j
                    sparsematrix(j,find(position==position(j)-vertice)) = -1;
                    neighbourtype(j,4) = find(position==position(j)-vertice);
                end
            else
                diagonal = diagonal -1; %if j is at x-edge diagonal is decreased
                xplanebot(end+1)=j;
            end
            %z-axis check in positive direction
            if position(j) <= depthnumber %if not at the max depth
                if ismember(position(j)+verticesq,position)==1 %is there something behind j
                    sparsematrix(j,find(position==position(j)+verticesq)) = -1;
                    neighbourtype(j,5) = find(position==position(j)+verticesq);
                end
            else
                diagonal = diagonal -1; %if j is at max depth diagonal is decreased
                zplanebot(end+1) =j;
            end
            %z-axis check in negative direction
            if position(j) > vertice^2 %if not at the front
                if ismember(position(j)-verticesq,position)==1 %is there something in front of j
                    sparsematrix(j,find(position==position(j)-verticesq)) = -1;
                    neighbourtype(j,6) = find(position==position(j)-verticesq);
                end
            else
                diagonal = diagonal -1; %if j is at front diagonal is decreased
                zplanetop(end+1)=j;
            end
            sparsematrix(j,j) = diagonal;
        end
        %making the matrix sparse after filling decreases computation time.
        sparsematrix = sparse(sparsematrix);
        %calculating the k value per voxel
        k = sparsematrix\ones(length(position),1);
        %generating a spaceholder matrix that will store series k values.
        keff = zeros(length(k),6);
        
        %% Second step with pressure gradient
        %for every k the neighbouring k's have their series k calculated
        %for the k_effective matrix.
        for j=1:length(k)
            valuek = k(j);
            neighbours = find(sparsematrix(j,:)==-1);
            if neighbours ~= 0
                for h=1:length(neighbours)
                    keff(j,h) = 2*k(neighbours(h))*valuek / (k(neighbours(h))+valuek);
                    keffmatrix(j,neighbours(h)) = -keff(j,h);
                end
                keffmatrix(j,j) = sum(keff(j,:));
            else %this should never be called but is a safety measure nonetheless.
                sparsematrix(j,j)=0;
            end
        end
        %Here the boundary conditions are applied, depending on whether we
        %work on x,y or z the diagonal elements of the matrix are
        %increased at the border voxels.
        %one boundary is set to pressure(1), the other to pressure(2) and
        %no-flow conditions are removed.
        if xyz==1
            if isempty(xplanebot)==0 && isempty(xplanetop)==0
                keffmatrixX=keffmatrix;
                boundary=zeros(length(position),1);
                for q=1:length(xplanetop)
                    if neighbourtype(xplanetop(q),4)~= 0
                        keffmatrixX(xplanetop(q),xplanetop(q)) = keffmatrixX(xplanetop(q),xplanetop(q)) - keffmatrixX(xplanetop(q),neighbourtype(xplanetop(q),4));
                        boundary(xplanetop(q),1) = pressure(1)*-keffmatrixX(xplanetop(q),neighbourtype(xplanetop(q),4));
                    end
                end
                for q=1:length(xplanebot)
                    if neighbourtype(xplanebot(q),3)~= 0
                        keffmatrixX(xplanebot(q),xplanebot(q)) = keffmatrixX(xplanebot(q),xplanebot(q)) - keffmatrixX(xplanebot(q),neighbourtype(xplanebot(q),3));
                        boundary(xplanebot(q),1) = pressure(2)*-keffmatrixX(xplanebot(q),neighbourtype(xplanebot(q),3));
                    end
                end
                keffmatrixX = sparse(keffmatrixX);
                px = keffmatrixX\boundary;
                velocity = zeros(length(xplanetop),1);
                for q=1:length(xplanetop)
                    if neighbourtype(xplanetop(q),4)~=0
                        velocity(q,1) = (pressure(1)-px(xplanetop(q)))*-keffmatrixX(xplanetop(q),neighbourtype(xplanetop(q),4));    
                    end
                end
                permeabilityX = sum(velocity)/((pressure(1)-pressure(2))*vertice);
                
            else
                permeabilityX =0;
            end
        end
        if xyz==2
            if isempty(yplanebot)==0 && isempty(yplanetop)==0
                keffmatrixY=keffmatrix;
                boundary=zeros(length(position),1);
                for q=1:length(yplanetop)
                    if neighbourtype(yplanetop(q),1)~= 0
                        keffmatrixY(yplanetop(q),yplanetop(q)) = keffmatrixY(yplanetop(q),yplanetop(q)) - keffmatrixY(yplanetop(q),neighbourtype(yplanetop(q),1));
                        boundary(yplanetop(q),1) = pressure(1)*-keffmatrixY(yplanetop(q),neighbourtype(yplanetop(q),1));
                    end
                    
                end
                for q=1:length(yplanebot)
                    if neighbourtype(yplanebot(q),2)~= 0
                        keffmatrixY(yplanebot(q),yplanebot(q)) = keffmatrixY(yplanebot(q),yplanebot(q)) - keffmatrixY(yplanebot(q),neighbourtype(yplanebot(q),2));
                        boundary(yplanebot(q),1) = pressure(2)*-keffmatrixY(yplanebot(q),neighbourtype(yplanebot(q),2));
                    end
                end
                keffmatrixY = sparse(keffmatrixY);
                py = keffmatrixY\boundary;
                velocity = zeros(length(yplanetop),1);
                for q=1:length(yplanetop)
                    if neighbourtype(yplanetop(q),1)~=0
                        velocity(q,1) = (pressure(1)-py(yplanetop(q)))*-keffmatrixY(yplanetop(q),neighbourtype(yplanetop(q),1));    
                    end
                end
                permeabilityY = sum(velocity)/((pressure(1)-pressure(2))*vertice);
                
            else
                permeabilityY =0;
            end
        end
        if xyz==3
            if isempty(zplanebot)==0 && isempty(zplanetop)==0
                keffmatrixZ=keffmatrix;
                boundary=zeros(length(position),1);
                for q=1:length(zplanetop)  
                    if neighbourtype(zplanetop(q),5)~= 0
                        keffmatrixZ(zplanetop(q),zplanetop(q)) = keffmatrixZ(zplanetop(q),zplanetop(q)) - keffmatrixZ(zplanetop(q),neighbourtype(zplanetop(q),5));
                        boundary(zplanetop(q),1) = pressure(1)*-keffmatrixZ(zplanetop(q),neighbourtype(zplanetop(q),5));
                    end
                    
                end
                for q=1:length(zplanebot)
                    if neighbourtype(zplanebot(q),6)~= 0
                        keffmatrixZ(zplanebot(q),zplanebot(q)) = keffmatrixZ(zplanebot(q),zplanebot(q)) - keffmatrixZ(zplanebot(q),neighbourtype(zplanebot(q),6));
                        boundary(zplanebot(q),1) = pressure(2)*-keffmatrixZ(zplanebot(q),neighbourtype(zplanebot(q),6));
                    end
                end
                keffmatrixZ = sparse(keffmatrixZ);
                pz = keffmatrixZ\boundary;                
                velocity = zeros(length(zplanetop),1);
                for q=1:length(zplanetop)
                    if neighbourtype(zplanetop(q),5)~=0
                        velocity(q,1) = (pressure(1)-pz(zplanetop(q)))*-keffmatrixZ(zplanetop(q),neighbourtype(zplanetop(q),5));    
                    end
                end
                permeabilityZ = sum(velocity)/((pressure(1)-pressure(2))*vertice);
                
            else
                permeabilityZ = 0;
            end
        end
    end
    perm(v,1) = (permeabilityX + permeabilityY + permeabilityZ)/3;
end
toc

%if the dataset is too large, use '-v7.3'
save('data\porosities.mat','percentages')
save('data\permeabilities.mat','perm')


%% Plotting
%uncomment to plot, though the plots only work for the last orientation
%that was calculated as this has the proper geometry.

% figure(1)
% fig = double(geometry);
% fig(position(:)) = px(:);
% slice(fig,1:vertice,1:vertice,1:vertice);
% 
% figure(2)
% fig = double(geometry);
% fig(position(:)) = py(:);
% slice(fig,1:vertice,1:vertice,1:vertice);

% figure(3)
% fig = double(geometry);
% fig(position(:)) = pz(:);
% slice(fig,1:vertice,1:vertice,1:vertice);

% figure(4)
% slice(double(geometry),1:vertice,1:vertice,1:vertice);

%% Questions?
% s.w.penninga@student.tue.nl