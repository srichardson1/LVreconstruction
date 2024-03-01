function Disp = DisplacementHex(nodeCurx, nodeCury, nodeCurz, ...
                           nodeRefx, nodeRefy, nodeRefz)
                       
Disp(1,:)=nodeCurx(:)'-nodeRefx(:)';
Disp(2,:)=nodeCury(:)'-nodeRefy(:)';
Disp(3,:)=nodeCurz(:)'-nodeRefz(:)';
