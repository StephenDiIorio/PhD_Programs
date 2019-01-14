import sys
sys.path.insert(0, '../../Utilities/')

from Particle import Particle
# from ParticleCloud import ParticleCloud
from Field import Field
import IOHelper as ioh
import OsirisConversions as oc
import numpy as np


global ndims, x_dim, y_dim, z_dim
ndims = 3
x_dim = 0
y_dim = 1
z_dim = 2
global nx, ny, nz, step, n0


def init_pos(npart):
    # init to their actual physical location. the grid grid and step size
    # are taken into account when we do the field interpolation
    centerx = 0.
    centery = ny * step[y_dim] / 2.
    centerz = nz * step[z_dim] / 2.

    spready = oc.pos_si_to_osiris_n0(15e-6, n0)
    spreadz = oc.pos_si_to_osiris_n0(15e-6, n0)

    to_ret = np.zeros((npart, ndims), dtype=np.float64)
    for index, pos in enumerate(to_ret):
        pos[x_dim] = centerx
        pos[y_dim] = (2 * spready) * np.random.random_sample() + \
                     (centery - spready)
        pos[z_dim] = (2 * spreadz) * np.random.random_sample() + \
                     (centerz - spreadz)
    return to_ret


def init_mom(npart):
    # to_ret = np.zeros(ndims, dtype=np.float64)
    # to_ret[x_dim] = np.float64(0.8035966691084001)  # 1E-3)
    lower_mom_tot = np.float64(0.6410483775684676)
    upper_mom_tot = np.float64(0.8035966691084001)

    # Corresponds to 1eV
    lower_py = 0.0  # np.float64(-0.00242294)
    upper_py = 0.0  # np.float64(0.00242294)
    lower_pz = 0.0  # np.float64(-0.00242294)
    upper_pz = 0.0  # np.float64(0.00242294)

    to_ret = np.zeros((npart, ndims), dtype=np.float64)
    momenta = (upper_mom_tot - lower_mom_tot) * \
        np.random.random_sample(size=npart) + lower_mom_tot
    for index, mom in enumerate(to_ret):
        py = (upper_py - lower_py) * np.random.random_sample() + lower_py
        pz = (upper_pz - lower_pz) * np.random.random_sample() + lower_pz

        px = np.sqrt(np.power(momenta[index], 2, dtype=np.float64) -
                     np.power(py, 2, dtype=np.float64) -
                     np.power(pz, 2, dtype=np.float64))

        mom[x_dim] = px
        mom[y_dim] = py
        mom[z_dim] = pz
    return to_ret


def init_efield(nx, ny, nz, ndims, strength, rmin, rmax):
    # e_field_amp = np.float64(0.000246127)
    e_field_amp = oc.efield_si_to_osiris_n0(strength, n0)
    print(e_field_amp)
    to_ret = np.zeros((ndims, nx, ny, nz), dtype=np.float64)

    for i in np.ndindex(nx, nz):
        ix, iz = i[0], i[1]
        centerx = ix - (nx / 2)
        centerz = iz - (nz / 2)
        radius = np.sqrt(np.power(centerx, 2) + np.power(centerz, 2))
        if radius >= rmin and radius < rmax:
            to_ret[x_dim, ix, :, iz] = np.cos(np.arctan2(centerz, centerx))
            to_ret[z_dim, ix, :, iz] = np.sin(np.arctan2(centerz, centerx))
    return to_ret * e_field_amp


def init_efield_fall(nx, ny, nz, ndims, strength, rmin, rmax):
    # e_field_amp = np.float64(0.000246127)
    e_field_amp = oc.efield_si_to_osiris_n0(strength, n0)
    print(e_field_amp)
    to_ret = np.zeros((ndims, nx, ny, nz), dtype=np.float64)

    for i in np.ndindex(nx, nz):
        ix, iz = i[0], i[1]
        centerx = ix - (nx / 2)
        centerz = iz - (nz / 2)
        radius = np.sqrt(np.power(centerx, 2) + np.power(centerz, 2))
        if radius >= rmin:
            to_ret[x_dim, ix, :, iz] = np.cos(np.arctan2(centerz, centerx)) * (rmin / radius)
            to_ret[z_dim, ix, :, iz] = np.sin(np.arctan2(centerz, centerx)) * (rmin / radius)
    return to_ret * e_field_amp


def init_bfield(nx, ny, nz, ndims, strength, rmin):
    to_ret = np.zeros((ndims, nx, ny, nz), dtype=np.float64)
    # to_ret[z_dim, :, :, :] = np.float64(1)
    return to_ret


def init_fields_zero(nx, ny, nz, ndims, strength, rmin):
    to_ret = np.zeros((ndims, nx, ny, nz), dtype=np.float64)
    return to_ret


#################################
# Main Functions
##################################


def single_run():
    global nx, ny, nz, step, n0

    np.random.seed(42)
    n0 = 0.17863390738e26  # density of simulated plasma [m^-3]

    # t = np.float64(0.0)
    dt = np.float64(0.1)  # 6)
    # tmax = np.float64(24 * np.pi)

    npart = np.uint16(10)

    pos_list = init_pos(npart)
    mom_list = init_mom(npart)

    nx = np.uint16(635)
    ny = np.uint16(5)
    nz = np.uint16(635)

    step = np.array([0.1, 0.1, 0.1], dtype=np.float64)
    inv_step = np.float64(1.0) / step

    # particles = ParticleCloud(npart, init_pos, init_mom)

    EField = Field(nx, ny, nz, init_efield)
    BField = Field(nx, ny, nz, init_bfield)

    # ioh.print_field_to_file(EField, save_x=True, save_y=False, save_z=True)

    # count = 0
    # data = np.empty((npart, int(tmax / dt) + 1, ndims), dtype=np.float64)
    data = [[] for i in range(npart)]

    # the max distance in x, using physical units, the particles should travel
    # this includes the distance through the field.
    max_x_dist = 143002 + (nx * step[0] / 2)
    # max_x_dist = 1000 + (nx * step[0] / 2)

    for i in range(npart):
        part = Particle(pos_list[i], mom_list[i], id_num=i)

        while not part.static:
            if part.in_domain:
                local_e_force, in_domain = \
                    EField.Weighted_Force(part.position,
                                          nx, ny, nz, inv_step)
                local_b_force, in_domain = \
                    BField.Weighted_Force(part.position,
                                          nx, ny, nz, inv_step)
                part.push_particle(local_e_force, local_b_force, dt)
                part.in_domain = in_domain
            else:
                part.single_particle_push(max_x_dist)
                part.static = True

            data[i].append(np.copy(part.position))

    # while particles.get_push_status() and t <= tmax:
    #     particles.push_cloud(nx, ny, nz, EField, BField, dt,
    #                          max_x_dist, inv_step)
    #     particles.update_push_status()

    #     for i in xrange(npart):
    #         data[i, count] = particles.get_particle_positions()[i]

    #     count += 1
    #     t += dt

    # data = np.delete(data, range(count, int(tmax / dt) + 1),
    #                  axis=1)  # remove unused data positions

    # gather data for the particles final position at a constant x value
    final_data = []
    for i in range(npart):
        final_data.append(data[i][-1])
    final_data = np.array(final_data)
    # final_data = data[:, -1, :]  # get the last position of the particles
    final_data = np.delete(final_data, 0, axis=1)  # remove x coord

    y = np.delete(final_data, 1, axis=1).flatten()
    z = np.delete(final_data, 0, axis=1).flatten()

    ioh.plot_trajectories(data)
    ioh.plot_heatmap(y, z)

    ioh.show_plots()


def multiple_sample():
    global nx, ny, nz, step, n0

    n0 = 0.17863390738e26  # density of simulated plasma [m^-3]

    # t = np.float64(0.0)
    dt = np.float64(0.1)  # 6)
    # tmax = np.float64(24 * np.pi)

    step = np.array([0.1, 0.1, 0.1], dtype=np.float64)
    inv_step = np.float64(1.0) / step

    npart = np.uint16(1000)

    shell_width = 0.1e-6

    radii = np.linspace(1e-6, 32e-6, 100)
    time = np.linspace(0., 100., 100)
    e0 = 1e9
    r0 = 1e-6
    strength = (e0 * r0) / np.array(radii)
    particle_count = []
    capture_range = 5.
    # data_range = [[7.5, 17.5], [-85, 115]]

    for i, r in enumerate(radii):
        np.random.seed(42)

        rmax = oc.pos_si_to_osiris_n0(r + shell_width, n0) * inv_step[x_dim]
        r = oc.pos_si_to_osiris_n0(r, n0) * inv_step[x_dim]
        print(r, rmax)
        # ensure that there is at least one cell difference between the two
        r = np.uint16(np.floor(r))
        rmax = np.uint16(np.ceil(rmax))
        print(r, rmax)

        nx = rmax
        nx = ((nx + 1) * 2)  # add one to account for rounding down
        nx = np.uint16(nx)

        ny = oc.pos_si_to_osiris_n0(15e-6, n0)
        ny = (((ny + 1) * 2) * inv_step[y_dim])
        ny = np.uint16(ny)

        nz = rmax
        nz = ((nz + 1) * 2)  # add one to account for rounding down
        nz = np.uint16(nz)

        data_range = [[10, 20], [((nz * step[z_dim]) / 2.) - capture_range, ((nz * step[z_dim]) / 2.) + capture_range]]

        pos_list = init_pos(npart)  # reset positions
        mom_list = init_mom(npart)  # reset momenta

        print(nx, ny, nz)

        EField = Field(nx, ny, nz, strength[i], r, rmax, init_efield)
        # BField = Field(nx, ny, nz, strength[i], r, init_bfield)

        # ioh.print_field_to_file(EField, save_x=True, save_y=False, save_z=True)

        data = [[] for j in range(npart)]

        # the max distance in x the particles should travel
        # this includes the distance through the field.
        max_x_dist = oc.pos_si_to_osiris_n0(0.18, n0) + (nx * step[x_dim] / 2)
        # max_x_dist = 1000 + (nx * step[x_dim] / 2)

        for j in range(npart):
            part = Particle(pos_list[j], mom_list[j], id_num=j)

            while not part.static:
                if part.in_domain:
                    local_e_force, in_domain = \
                        EField.Weighted_Force(part.position,
                                              nx, ny, nz, inv_step)

                    part.push_particle_no_b(local_e_force, dt)
                    part.in_domain = in_domain
                else:
                    part.single_particle_push(max_x_dist)
                    part.static = True

                data[j].append(np.copy(part.position))

        # gather data for the particles final position at a constant x value
        final_data = []
        for j in range(npart):
            final_data.append(data[j][-1])
        final_data = np.array(final_data)
        final_data = np.delete(final_data, 0, axis=1)  # remove x coord
        y = final_data[:, 0]
        z = final_data[:, 1]

        filename = str(i) + "_r=" + str(radii[i]) + "_ef=" + str(strength[i])
        ioh.save_data_array(filename, final_data)
        ioh.plot_trajectories(data, filename + '_traj')
        h, h_r = ioh.plot_heatmap(y, z, data_range, filename + '_heat')

        particle_count.append(np.sum(h_r))
        print(particle_count)

        # ioh.show_plots()
        ioh.close_figures()

    b_y = (30 - 0) * np.random.random_sample(npart) + 0
    b_z = (300 - (-300)) * np.random.random_sample(npart) + (-300)
    b_heatmap, yedges, zedges = np.histogram2d(b_y, b_z, bins=100)
    b_heatmap_r, yedges_r, zedges_r = np.histogram2d(b_y, b_z, range=data_range, bins=100)
    base_count = np.sum(b_heatmap_r)

    ioh.plot_rel_diff(radii, time, particle_count, base_count, str(i) + '_rel_data')
    ioh.plot_count_diff(radii, time, particle_count, str(i) + '_count_data')
    ioh.show_plots()

#################################


if __name__ == "__main__":
    # single_run()

    multiple_sample()
