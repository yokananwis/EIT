function [X,Y] = coord_for_elec(height,no,r)
% 
% Laskee karteesiset koordinaatit 2d elektrodeille, kun
% elektrodien koko, lukumaara ja tankin sade tiedetaan.
%
% L. Heikkinen 6.5.1999
%

elecA = height*no;
p = 2*pi*r;
otherA = p - elecA;
wall = otherA/no;

elecangle = (height*pi)/(pi*r); 
ele = elecangle/3;
wallangle = (wall*pi)/(pi*r);

k = wallangle*ones(1,no*3+no);

k(1:4:no*3+no) = ele;
k(2:4:no*3+no) = ele;
k(3:4:no*3+no) = ele;

k = diag(k);

m = tril(ones(no*3+no));

l = m*k;
angle = sum(l,2);

[x,y] = pol2cart(angle,r);
s = size(x,1);
X = [x(s);x(1:s-1)];
s = size(y,1);
Y = [y(s);y(1:s-1)];
