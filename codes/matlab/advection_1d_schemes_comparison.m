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
%  Code: 1D advection transport equation solved using upwind scheme       %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% prepare video
v = VideoWriter('advection_1d_comparison.mp4', 'MPEG-4');
open(v);

% one-dimensional advection by first order upwind
n=41;
dt=0.05;
u=1;
h=4/(n-1);
x=1:n;

fanalytical=zeros(n,1);
yupwind=zeros(n,1);
fupwind=zeros(n,1);
yleap=zeros(n,1);
yleapold=zeros(n,1);
fleap=zeros(n,1);
ylax=zeros(n,1);
flax=zeros(n,1);
ymac=zeros(n,1);
fmac=zeros(n,1);

fanalytical(1)=1.0;
fupwind(1)=1.0;
fleap(1)=1.0;
flax(1)=1.0;
fmac(1)=1.0;
nstep=56;
Co = u*dt/h;

t=0;
for m=1:nstep
    
    fprintf('Step %d: - Courant number: %f - Time: %f\n', m, Co, t);
    
    hold off;
    subplot(221);
    plot(x,fanalytical, x, fupwind);
    axis([1, n, -0.5, 1.5]);
    title('Upwind');
    
    hold off;
    subplot(222);
    plot(x,fanalytical, x, fleap);
    axis([1, n, -0.5, 1.5]);
    title('Leap Frog');
    
    hold off;
    subplot(223);
    plot(x,fanalytical, x, flax);
    axis([1, n, -0.5, 1.5]);
    title('Lax-Wendroff');
    
    hold off;
    subplot(224);
    plot(x,fanalytical, x, fmac);
    axis([1, n, -0.5, 1.5]);
    title('MacCormack');
    
    pause(0.025);
    
    frame = getframe(gcf);
    writeVideo(v,frame);

    yupwind=fupwind;
    yleapold=yleap;
    yleap=fleap;
    ylax=flax;
    ymac=fmac;

    for i=2:n-1
        
        fanalytical(i) = (sign(u*t-h*(i-1))+1)/2;
        
        fupwind(i)=yupwind(i)-(u*dt/h)*(yupwind(i)-yupwind(i-1));
        
        fleap(i)=yleapold(i)-(u*dt/h)*(yleap(i+1)-yleap(i-1));
        
        flax(i)=ylax(i)-(u*dt/2/h)*(ylax(i+1)-ylax(i-1)) + ...
                u*u*dt*dt/2/h/h*(ylax(i+1)-2*ylax(i)+ylax(i-1));
        
        ystarc = ymac(i)-(u*dt/h)*(ymac(i+1)-ymac(i));
        ystarb = ymac(i-1)-(u*dt/h)*(ymac(i)-ymac(i-1));
        fmac(i)=0.5*(ymac(i)+ystarc-(u*dt/h)*(ystarc-ystarb));
        
    end
    
    t = t+dt;
end

close(v);

