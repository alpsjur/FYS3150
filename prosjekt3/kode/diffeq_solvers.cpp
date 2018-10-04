

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
