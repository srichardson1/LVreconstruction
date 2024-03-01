external_two_left=zeros(65537,1);
external_left=zeros(65537,1);
external_middle=zeros(65537,1);
external_right=zeros(65537,1);
external_two_right=zeros(65537,1);

for j=1:65537
    external_two_right(j,1)=external(j,1);
    external_right(j,1)=external(j+65537+1211,1);
    external_middle(j,1)=external(j+65537+2422,1);
    external_left(j,1)=external(j+65537+3633,1);
    external_two_left(j,1)=external(j+65537+4844,1);
end

figure
hold on
plot(time(1:65537),external_two_right)
plot(time(1:65537),aortic_new)

figure
hold on
plot(time(1:65537),external_right)
plot(time(1:65537),aortic_new)

figure
hold on
plot(time(1:65537),external_middle)
plot(time(1:65537),aortic_new)

figure
hold on
plot(time(1:65537),external_left)
plot(time(1:65537),aortic_new)


figure
hold on
plot(time(1:65537),external_two_left)
plot(time(1:65537),aortic_new)