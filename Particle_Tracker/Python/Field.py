import numpy as np
from FrozenClass import FrozenClass

ndims = np.uint8(3)


class Field(FrozenClass):
    def __init__(self, nx, ny, nz, strength, rmin, rmax, f):
        self.field = f(nx, ny, nz, ndims, strength, rmin, rmax)
        return

    def set_field(self, f):
        msg = ("incompatible dimensions for field\n"
               "(dimension must be %r)")
        if f.shape != self.field.shape:
            raise ValueError(msg % (self.field.shape,))

        self.field = f
        return

    def get_field(self):
        return self.field

    def Weighted_Force(self, pos, nx, ny, nz, inv_step_size):
        # linear interpolation in 3d

        # variables to check whether any of the particle still lies
        # within the field domain
        domain_count = 0
        domain_status = True

        # charge to mass ratio - normalized to |e|/m_e
        qoverm = np.float64(-1.0)

        local_force = np.zeros(ndims, dtype=np.float64)

        volumes = np.zeros((2, 2, 2), dtype=np.float64)

        weight = np.copy(pos)
        weight *= inv_step_size
        ii = weight.astype(np.int64)
        weight -= ii.astype(np.float64)

        for index, in np.ndenumerate(volumes):
            volumes[index] = np.absolute(((np.float64(1.0) - index[0]) - weight[0]) * ((np.float64(1.0) - index[1]) - weight[1]) * ((np.float64(1.0) - index[2]) - weight[2]), dtype=np.float64)

        for index, in np.ndenumerate(volumes):
            i0 = ii[0] + index[0]
            i1 = ii[1] + index[1]
            i2 = ii[2] + index[2]
            if ((i0 < nx) and (i1 < ny) and (i2 < nz) and
                    (i0 >= 0) and (i1 >= 0) and (i2 >= 0)):
                local_force += volumes[index] * self.field[:, i0, i1, i2]
            else:
                domain_count += 1
                # print "WARNING: Evaluating field outside specified domain. (It treats this as if there is no field at this position.)"

        if domain_count == volumes.size:
            domain_status = False

        return local_force * qoverm, domain_status
