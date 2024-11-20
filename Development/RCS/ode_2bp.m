function dy = ode_2bp(~, y,  mu)

r = y(1:3);
v = y(4:6);

dy = [v; (-mu/(norm(r)^3)) * r];

end

