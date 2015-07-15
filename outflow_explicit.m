function prim_cc_go = outflow_explicit(qN,pback,CSA_flag)

if nargin == 2
    CSA_flag = 0;
end

if pback<0
    if CSA_flag ~=1
        prim_cc_go = qN;
    else
        % For C-D nozzle
        prim_cc_go = qN;
%         prim_cc_go = [ -0.247429426601287;
%                        1.997508077914113e+002;
%                        -1.978124262577476e+004];
    end
else
    if CSA_flag ~=1
        prim_cc_go = qN;
        prim_cc_go(3) = pback;
    else
        prim_cc_go = [ 0.009961807199782;
                      -6.082036704878391;
                       0.0];
%         prim_cc_go = 2*qN - qNm1;
%         prim_cc_go(3) = 0;
    end
end
end