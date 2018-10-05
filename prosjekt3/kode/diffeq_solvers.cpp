

void velocity_verlet(double dt, double *acc, double *pos, double *vel){
  double dt2 = dt/2.0;
  double 2dt2 = dt*dt/2.0;
  double acc_new;


  for(int dim = 0; dim < 3; ++dim){
    pos[dim] = pos[dim] + dt*vel[dim] + 2dt2*acc[dim];
    acc_new = calculate_acc(pos[dim]);
    vel[dim] = vel[dim] + dt2*(acc_new + acc[dim]);
    acc[dim] = acc_new;
  }
}

void velocityVerlet(double dt, Coordinate &acc, Coordinate &pos, Coordinate &vel){
  double dt = dt/2.0;
  double 2dt2 = dt*dt/2.0;
  Coordinate acc_new;

  pos = pos + dt*vel + 2dt2*acc;
  acc_new = calculate_acc(pos);
  vel = vel + dt2*(acc_new + acc);
  acc = acc_new; 
}

void forwardEuler(double dt, Coordinate &pos, Coordinate &vel){
  Coordinate acc = calculate_acc(pos);
  pos = pos + vel*dt;
  vel = vel + acc*dt;
}
