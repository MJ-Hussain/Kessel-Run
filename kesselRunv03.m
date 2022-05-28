clc;
clear;
%Time step
dt=0.1;

%Number of time steps
T=2000;

%mass of ship
ms=1;

%Reading black hole cluster position and mass data
file=load("cluster1.mat");
% mass data
mg=file.hM;
%position data
pg(:,1)=file.hX;
pg(:,2)=file.hY;

%Number of black holes
m_len=length(mg);

%Counter fot successfull run
sc=0;

%Number of histories to run
N=100;

%Running loop for histories
for j=1:N
  
%initial position of ship
% Random position on X-axis b/w -5 to 5
prx=unifrnd(-5,5);
%position vector
pr=[prx,-10];

%initial velocity magnitude of the ship
v_mag=unifrnd(2,5);
% Normal random direction of ship
v_dir=normrnd(pi/2,pi/4);
% Velocity vector
v(1)=v_mag*cos(v_dir);
v(2)=v_mag*sin(v_dir);

% Loop for time steps
  for t=2:T+1
    
    %Initializing force vector
    F=[0,0];
    
    %Net force from all black holes
    for g=1:m_len  
      %Position vector from black hole
      r=pr(t-1,:)-pg(g,:);
      %Positon vector magnitude
      r_norm=norm(r);
      %Calculating Force
      F= F+(-(ms*mg(g)/r_norm^3)*r);
    end
  
  %Acceleration
  a=F/ms;
  
  %Condition on acceleration 
  %Loop termination when acceleration exceeds 4 m/s^2
  if (norm(a)>=4)
    break
    
  else
    %Velocity at current time step
    v=v+dt*a;
    %New position from updated velocity
    pr(t,:)=pr(t-1,:)+dt*v;
  end
  
  %Termination when ship position exceeds 10 on Y-axis
  if (pr(t,2)>10)
    %successfull history added
    sc=sc+1;
    %Position of successfull history saved
    po(sc,1:t,:)=pr(:,:);
    %Time of successfull history saved
    time(sc)=t;
    %Successfull history initial position
    p_opt(sc)=prx;
    %Successfull history initial velocity
    v_opt(sc)=v_mag;
    break
    
  %Termination if ship position is > 10 or < -10 i.e. out of region of interest
  elseif ((pr(t,1)>10)||(pr(t,1)<-10))
    break
    
  end
    end
  end

%Maximum path time of ship 
[max_t ,ind_tmax]=max(time);
%Minimum path time of ship 
[min_t ,ind_tmin]=min(time);

%Longest path position vectors
long(:,:)=po(ind_tmax,:,:);
%Shortest path position vectors
short(:,:)=po(ind_tmin,:,:);
%Eliminating rows of zeros
ind = find(sum(short,2)==0) ;
short(ind,:) = [] ;

%values for color map
vec1=length(long(:,1));
%RGB triplet matrix
c1=zeros(vec1,3);
c1(:,1)=linspace(0,1,vec1);
c1(:,3)=linspace(1,0,vec1);

vec2=length(short(:,1));
%RGB triplet matrix
c2=zeros(vec2,3);
c2(:,1)=linspace(0,1,vec2);
c2(:,3)=linspace(1,0,vec2);

fprintf("optimum starting position of the ship is: %3.1f\n",p_opt(ind_tmin))
fprintf("Optimal velocity of the ship is: %3.1f\n",v_opt(ind_tmin))

%ploting the Ship path using a color map scatter plot
%Shortest path
scatter(short(:,1),short(:,2),[], c2, "filled");
hold on

%Longest path
scatter(long(:,1),long(:,2),[], c1, "filled");
hold on


%Plotting black holes positions
scatter(pg(:,1), pg(:,2));

ylim([-20 20]);
xlim([-20 20]);
hold off