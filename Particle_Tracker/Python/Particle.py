import numpy as np
from FrozenClass import FrozenClass

ndims = np.uint8(3)


class Particle(FrozenClass):
    def __init__(self,
                 pos=np.zeros(ndims, dtype=np.float64),
                 mom=np.zeros(ndims, dtype=np.float64),
                 in_domain=True,
                 static=False,
                 id_num=1):
        msg = ("incompatible dimensions for %s\n"
               "(dimension must be %d)")
        if pos.shape[-1] != ndims:
            raise ValueError(msg % ('position', ndims))
        if mom.shape[-1] != ndims:
            raise ValueError(msg % ('momentum', ndims))

        self.position = pos
        self.momentum = mom
        self.in_domain = in_domain
        self.static = static
        self.id = id_num
        return

    def set_id(self, num):
        self.id = num
        return

    def get_id(self):
        return self.id

    def set_in_domain(self, in_domain):
        self.in_domain = in_domain
        return

    def is_in_domain(self):
        return self.in_domain

    def set_static(self, static):
        self.static = static
        return

    def is_static(self):
        return self.static

    def set_pos(self, pos):
        msg = ("incompatible dimensions for position\n"
               "(dimension must be %d)" % (ndims,))
        if pos.shape[-1] != ndims:
            raise ValueError(msg)

        self.position = pos
        return

    def get_pos(self):
        return self.position

    def set_mom(self, mom):
        msg = ("incompatible dimensions for momentum\n"
               "(dimension must be %d)" % (ndims,))
        if mom.shape[-1] != ndims:
            raise ValueError(msg)

        self.momentum = mom
        return

    def get_mom(self):
        return self.momentum

    def single_particle_push(self, dist_x):
        # mom2 = np.dot(self.momentum, self.momentum)
        # gamma = np.float64(1.0) / np.sqrt(np.float64(1.0) + mom2)

        t = (dist_x - self.position[0]) / self.momentum[0]  # the x-comp

        self.position += self.momentum * t  # * gamma

        return

    def push_particle(self, e_force, b_force, dt):
        msg = ("incompatible dimensions for %s\n"
               "(dimension must be %d)")
        if e_force.shape[-1] != ndims:
            raise ValueError(msg % ('electric field force', ndims))
        if b_force.shape[-1] != ndims:
            raise ValueError(msg % ('magnetic field force', ndims))

        dt = np.float64(dt)

        self.__push_momentum(e_force, dt * np.float64(0.5))
        self.__lorentz(b_force, dt)
        self.__push_momentum(e_force, dt * np.float64(0.5))
        self.__push_pos(dt)
        return


    def push_particle_no_b(self, e_force, dt):
        msg = ("incompatible dimensions for %s\n"
               "(dimension must be %d)")
        if e_force.shape[-1] != ndims:
            raise ValueError(msg % ('electric field force', ndims))

        dt = np.float64(dt)

        self.__push_momentum(e_force, dt * np.float64(0.5))
        self.__push_momentum(e_force, dt * np.float64(0.5))
        self.__push_pos(dt)
        return


    def __push_momentum(self, local_e, dt):
        dt = np.float64(dt)

        self.momentum += local_e * dt
        return

    def __push_pos(self, dt):
        dt = np.float64(dt)

        mom2 = np.dot(self.momentum, self.momentum)
        gamma = np.float64(1.0) / np.sqrt(np.float64(1.0) + mom2)

        self.position += self.momentum * dt * gamma
        return

    def __lorentz(self, local_b, dt):
        dt = np.float64(dt)

        loc_b2 = np.dot(local_b, local_b)

        if loc_b2 != np.float64(0.0):
            mom2 = np.dot(self.momentum, self.momentum)
            gamma = np.float64(1.0) / np.sqrt(np.float64(1.0) + mom2)

            tt2 = np.float64(0.0)
            tt = local_b * dt * np.float64(0.5)
            tt2 += np.dot(tt, tt)

            ss = np.float64(2.0) * tt / (np.float64(1.0) + tt2)

            vperp = self.momentum * (np.float64(1.0) - local_b / loc_b2)
            vstar = vperp + np.cross(vperp, tt) * gamma

            self.momentum += np.cross(vstar, ss) * gamma
        return
