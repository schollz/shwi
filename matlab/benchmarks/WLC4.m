% function [F] = WLC(p,x,L)
%     F = 4.1 / p * ( ones(size(x)) ./ (4 * (1-x./L).*(1-x./L)) + [x./L] - 1/4);
% end
% %Estimating the Persistence Length of a Worm-Like Chain Molecule from
%Force-Extension Measurements
% Croquette
function [F] = WLC4(p,L,x)
%     F = 4.1 / p * ( ones(size(x)) ./ (4 * (1-x./L).*(1-x./L)) + [x./L] - 1/4);
   F = 4.1 / p * ( ones(size(x)) ./ (4 * (1-x./L).*(1-x./L)) + [x./L] - 1/4 + -0.5164228*[x./L].^2 -2.737418*[x./L].^3 + 16.07497*[x./L].^4 -38.87607*[x./L].^5+39.49944*[x./L].^6 -14.17718*[x./L].^7);
     F(x>L)=0;
end
