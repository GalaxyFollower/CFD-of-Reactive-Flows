% ----------------------------------------------------------------------- %
%     ╔═══╦═══╦═══╗   ╔═╦═══╗        ╔╗       ╔═══╦╗                      %
%     ║╔═╗║╔══╩╗╔╗║   ║╔╣╔═╗║       ╔╝╚╗      ║╔══╣║                      %
%     ║║─╚╣╚══╗║║║╠══╦╝╚╣╚═╝╠══╦══╦═╩╗╔╬╦╗╔╦══╣╚══╣║╔══╦╗╔╗╔╦══╗          %
%     ║║─╔╣╔══╝║║║║╔╗╠╗╔╣╔╗╔╣║═╣╔╗║╔═╣║╠╣╚╝║║═╣╔══╣║║╔╗║╚╝╚╝║══╣          %
%     ║╚═╝║║  ╔╝╚╝║╚╝║║║║║║╚╣║═╣╔╗║╚═╣╚╣╠╗╔╣║═╣║  ║╚╣╚╝╠╗╔╗╔╬══║          % 
%     ╚═══╩╝  ╚═══╩══╝╚╝╚╝╚═╩══╩╝╚╩══╩═╩╝╚╝╚══╩╝  ╚═╩══╝╚╝╚╝╚══╝          %
% ----------------------------------------------------------------------- %
%                                                                         %
%   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       %
%   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            %
%   Department of Chemistry, Materials and Chemical Engineering           %
%   Politecnico di Milano                                                 %
%   P.zza Leonardo da Vinci 32, 20133 Milano                              %
%                                                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
%   This file is part of CFDofReactiveFlows framework.                    %
%                                                                         %
%   License                                                               %
%                                                                         %
%   Copyright(C) 2018 Alberto Cuoci                                       %
%   CFDofReactiveFlows is free software: you can redistribute it and/or   %
%   modify it under the terms of the GNU General Public License as        %
%   published by the Free Software Foundation, either version 3 of the    %
%   License, or (at your option) any later version.                       %
%                                                                         %
%   CFDofReactiveFlows is distributed in the hope that it will be useful, %
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%   GNU General Public License for more details.                          %
%                                                                         %
%   You should have received a copy of the GNU General Public License     %
%   along with CFDofReactiveFlows.                                        %
%   If not, see <http://www.gnu.org/licenses/>.                           %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%  Code: 1D advection transport equation showing the concept of           %
%        artificial viscosity                                             %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% prepare video
v = VideoWriter('artificial_viscosity.mp4', 'MPEG-4');
open(v);

% data
n=161; 
nstep=250; 
length=4.0;
h=length/(n-1);
dt=0.25*h;
u=1;
gamma=0.5;
Co = u*dt/h;

yupwind=zeros(n,1);
fupwind=zeros(n,1);
fupwind(1)=1;
ylax=zeros(n,1);
flax=zeros(n,1);
flax(1)=1;
ylaxart=zeros(n,1);
flaxart=zeros(n,1);
flaxart(1)=1;

% Iterations
t =0;
for m=1:nstep
    
    fprintf('Step %d: - Courant number: %f\n - Time: %f', m, Co, t);
    
    subplot(3,1,1);
    hold off;
    plot(yupwind,'linewidt',2); 
    axis([1 n -0.5, 1.5]);
    hold on;
    plot([1,dt*(m-1)/h+1.5,dt*(m-1)/h+1.5,n],[1,1,0,0],'r','linewidt',2);
    title('Upwind');
    
    subplot(3,1,2);
    hold off;
    plot(ylax,'linewidt',2); 
    axis([1 n -0.5, 1.5]);
    hold on;
    plot([1,dt*(m-1)/h+1.5,dt*(m-1)/h+1.5,n],[1,1,0,0],'r','linewidt',2);
    title('Lax-Wendroff');
    
    subplot(3,1,3);
    hold off;
    plot(ylaxart,'linewidt',2); 
    axis([1 n -0.5, 1.5]);
    hold on;
    plot([1,dt*(m-1)/h+1.5,dt*(m-1)/h+1.5,n],[1,1,0,0],'r','linewidt',2);
    title('Lax-Wendroff with D=0.5');
    
    pause(0.0125);
    frame = getframe(gcf);
    writeVideo(v,frame);

    yupwind=fupwind; 
    ylax=flax; 
    ylaxart=flaxart; 
    t=t+dt;

    
    for i=2:n-1
        
        fupwind(i) = yupwind(i)-(u*dt/h)*(yupwind(i)-yupwind(i-1));
                    
        flax(i) = ylax(i)-(0.5*u*dt/h)*(ylax(i+1)-ylax(i-1))+...
                          (0.5*(u*dt/h)^2)*(ylax(i+1)-2.0*ylax(i)+ylax(i-1));
                                    
        flaxart(i) = ylaxart(i)-(0.5*u*dt/h)*(ylaxart(i+1)-ylaxart(i-1))+...
                     (0.5*(u*dt/h)^2)*(ylaxart(i+1)-2.0*ylaxart(i)+ylaxart(i-1))+...
                     gamma*(dt/h)*( abs(ylaxart(i+1)-ylaxart(i))*(ylaxart(i+1)-ylaxart(i))-...
                                    abs(ylaxart(i)-ylaxart(i-1))*(ylaxart(i)-ylaxart(i-1))  );
    end
    
end

close(v);