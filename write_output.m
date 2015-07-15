function write_output(n, x_cc, area_cc, prim_cc, CSA_flag)

if nargin == 4
    CSA_flag = 0;
end

N = size(prim_cc,2);
if CSA_flag ~= 1
    fp2 = fopen('./nozzle.dat','a');
else
    fp2 = fopen('./nozzle_local.dat','a');
end
fprintf(fp2, 'zone T="n=%d"\n',n);
fprintf(fp2, 'I= %d\n',N);
fprintf(fp2, 'DATAPACKING=POINT\n');
for i=1:N
    fprintf(fp2,'%e %e %e %e %e\n', x_cc(i), area_cc(i), ...
        prim_cc(1,i), prim_cc(2,i), prim_cc(3,i));
end
fclose(fp2);
end