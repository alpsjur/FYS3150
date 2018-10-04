

class Planet{
  // constant/variable declarations
  const double m_mass;
  double * m_init_pos, m_init_vel [3];

  // function declarations
public:
  Planet(double mass, double * init_pos, double * init_vel)
    : m_mass(mass), m_init_pos(init_pos), m_init_vel(init_vel)
  {
    // initialising mass, initial position and initial velocity of planet
  }
}
