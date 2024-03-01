
PP=[36.0609  -49.6045   16.7145
   37.1406  -50.1238   15.4619
   42.9907  -56.5533   10.4725
   48.6689  -62.9674    5.2882
   53.3860  -71.1563    1.2469
   56.6470  -78.2393   -5.6173
   61.9193  -85.6682  -10.0006
   65.8621  -93.9632  -14.7098];
XX=PP(:,1);YY=PP(:,2);ZZ=PP(:,3);
figure(1);
scatter3(XX,YY,ZZ,'o');
xlabel('x');
ylabel('y');
zlabel('z');

%NormalizationVector
LAVec = NormalizationVec(PP(1,:) - PP(:,8));

SAXVec = NormalizationVec(DataSegSA(1,1).endo_cReal(:,1)'-LVUpperCenter);
    SAXVec_t = NormalizationVec(DataSegSA(1,1).endo_cReal(:,10)'-LVUpperCenter);
    LAVec_t = cross(SAXVec,SAXVec_t);
