function [Element,Nodelist]=MakeElement(H,bg,E);

%MakeElement Computes the Element data structure for a 2D mesh
% Function [Element,Nodelist]=MakeElement(H,bg,E); 
% computes the Element data for a mesh. 
% Element is a structure including the topology (H) of each element,
% faces of each element, adjacent element of each face and information
% whether the face is internal or boundary face.
% bg are the boundary node indices and Nodelist includes all the faces.   
% E includes the electrode information.
%
% INPUT
%
% H = topology matrix
% bg = indices of the boundary nodes
% E = elements under the electrodes having two nodes on the boundary
%
% OUTPUT
%
% Element = element data structure
% Nodelist = list of all the faces

% M. Vauhkonen 10.10.1998. Modified for EIDORS 11.8.1999 by 
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


Ind=[1 2;2 3;3 1];
ellist=[];
Nodelist=[];
[rH,cH]=size(H);
 for ii=1:rH
  Element(ii).Topology=[H(ii,:)]; %Topology tells which of the nodes 
                                  % are included in the ii'th element.
  for jj=1:3
   ind1=H(ii,Ind(jj,1));
   ind2=H(ii,Ind(jj,2));
   Nodelist=[Nodelist;[ind1,ind2]];
   [I1,I2]=find(ind1==H);
   I=[];
    for kk=1:size(I1,1)
      if find(H(I1(kk),:)==ind2)
       	I=[I,I1(kk)];
      end
    end
   I=sort(I);
   Element(ii).Face{jj,1}=[ind1,ind2]; % Indices of the nodes of the ii'th 
                                       % element's jj'th face.  
    if find(ind1==bg) & find(ind2==bg)
       [fE1,FE2]=find(ii==E);          % Face under the electrode?
       if ~isempty(fE1)
        Element(ii).Face{jj,3}=fE1;    % Insert electrode index.
       else
        Element(ii).Face{jj,3}=0;     
       end   
      Element(ii).Face{jj,2}=0;        % Adjacent element index. If zero, 
                                       % the face is on the boundary. 
    else   
      Element(ii).Face{jj,2}=I(I~=ii); % Adjacent element index.      
      Element(ii).Face{jj,3}=0;        % Not under any electrode.
    end
  end  	
 end 











