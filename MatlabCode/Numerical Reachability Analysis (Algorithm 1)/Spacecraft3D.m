classdef Spacecraft3D    
% -------------------------------------------------------------------------
%  This class descripbe the physical properties of a solar sail
%
% [inputs]: - [Array] : 7 by 1 vector composed of [L l cb M m rhos rhod Cr]
%
% [outputs]: - [I_total] : Total moment of inertia of the sailcraft
% -------------------------------------------------------------------------

   properties
     L    % Length of the solar sail
     l    % width of the solar sail 
     r    % Radius of the spacecraft (for cannonball model)
     cb   % Control boom
     M    % Mass of the bus
     m    % Mass of the coner supports
     rhos % rate of specular reflection
     rhod % rate of diffusion
     Cr   % reflectivity coefficient
     side % Side of the cubesat
     CD   % Drag coeficient
   end
   methods
       
       function Current_Spacraft = Spacecraft3D(Array)
           if length(Array) == 11
               Current_Spacraft.L         = Array(1);
               Current_Spacraft.l         = Array(2);
               Current_Spacraft.r         = Array(3);
               Current_Spacraft.cb        = Array(4);
               Current_Spacraft.M         = Array(5);
               Current_Spacraft.m         = Array(6);
               Current_Spacraft.rhos      = Array(7);
               Current_Spacraft.rhod      = Array(8);
               Current_Spacraft.side      = Array(9);
               Current_Spacraft.CD        = Array(10);
               Current_Spacraft.Cr        = Array(11);
           end
           
       end
       
      function I_total = Moment_Of_Inertia_Calculator(obj)
         mass = obj.M+obj.m;
         I_Flatplate = [1/12*mass*obj.L^2; 1/12*mass*obj.l^2; 1/12*mass*(obj.L^2+obj.l^2)];
         I_cube = [mass*obj.side^2/6; mass*obj.side^2/6; mass*(obj.side^2/6 + obj.side^2/2)];
         I_total = I_Flatplate + I_cube;
      end
      
   end
end















