import numpy as np
from FrozenClass import FrozenClass
from Particle import Particle


class ParticleCloud(FrozenClass):
    def __init__(self, num_part, pos_f, mom_f, push_status=True):
        self.particles = np.empty(num_part, dtype=object)

        pos_list = pos_f(num_part)
        mom_list = mom_f(num_part)
        for i in range(num_part):
            self.particles[i] = Particle(pos=pos_list[i],
                                         mom=mom_list[i],
                                         id_num=i)
        self.still_push = push_status
        return

    def set_push_status(self, status):
        self.still_push = status
        return

    def get_push_status(self):
        return self.still_push

    def get_particles(self):
        return self.particles

    def get_particle_positions(self):
        to_ret = np.empty_like(self.particles)

        for index, part in enumerate(self.particles):
            to_ret[index] = part.position

        return to_ret

    def get_particle_momenta(self):
        to_ret = np.empty_like(self.particles)

        for index, part in enumerate(self.particles):
            to_ret[index] = part.momentum

        return to_ret

    def update_push_status(self):
        still_pushing = False
        for index, part in enumerate(self.particles):
            if not part.static:
                still_pushing = True
                break
        self.still_push = still_pushing
        return still_pushing

    def push_cloud(self, nx, ny, nz, e_field, b_field,
                   dt, dist_x, inv_step_size):
        for index, part in enumerate(self.particles):
            if part.in_domain:
                # still in domain so do normal push
                local_e_force, in_domain = \
                    e_field.Weighted_Force(part.position,
                                           nx, ny, nz, inv_step_size)
                local_b_force, in_domain = \
                    b_field.Weighted_Force(part.position,
                                           nx, ny, nz, inv_step_size)

                part.push_particle(local_e_force, local_b_force, dt)
                part.in_domain = in_domain
            elif not part.static:
                # outside of domain, so if the particle is not done moving,
                # complete the final motion in a single time step (should
                # just be a linear trajectory)
                part.single_particle_push(dist_x)
                part.static = True
            else:
                # otherwise, the particle is done moving, so we don't need
                # bother looking at it anymore
                pass
        return
