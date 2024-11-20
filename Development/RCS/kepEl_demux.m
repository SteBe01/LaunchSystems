function [a, e, i, OM, om, f] = kepEl_demux(kepEl)
% It exports keplerian elements from vector
%
% PROTOTYPE:
% [a, e, i, OM, om, f] = kepEl_demux(kepEl)
%
% INPUT :
%    kepEI   [1x6]   = vector of keplerian elements
%
%
%    kepEI(1)        = a                              [km]
%    kepEI(2)        = e                              [-]
%    kepEI(3)        = i                              [rad]
%    kepEI(4)        = OM                             [rad]
%    kepEI(5)        = om                             [rad]
%    kepEI(6)        = theta                          [rad]ional parameter of the primary [L^3/T^2]
%
% CONTRIBUTORS:
% Camilla Airoldi, Jesus Braza Rodriguez, Elena De Marco, Riccardo Monti
%
% VERSIONS
% 2023-12-20 Last version  


a  = kepEl(1);
e  = kepEl(2);
i  = kepEl(3);
OM = kepEl(4);
om = kepEl(5);
f  = kepEl(6);

end