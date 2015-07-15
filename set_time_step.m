function dt = set_time_step(dx, prim_cc,CFL)
g=1.4;
a = sqrt(g*prim_cc(3,:)./prim_cc(1,:));

dt_temp = CFL*dx./(abs( a + prim_cc(2,:) ));

dt = min( min(dt_temp), 10 )*ones(size(dx));

end