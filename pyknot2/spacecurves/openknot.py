'''
OpenKnot
========

A class for working with the topology of open curves.
'''

from __future__ import print_function
import numpy as n
import signal
import copy

from pyknot2.spacecurves.spacecurve import SpaceCurve
from pyknot2.spacecurves.knot import Knot
from pyknot2.spacecurves.rotation import get_rotation_angles, rotate_to_top
from collections import Counter

from pyknot2.invariants import alexander, generalised_alexander, vassiliev_degree_2
from pyknot2.representations.gausscode import GaussCode

class Timeout:
  """Timeout class using ALARM signal"""
  class Timeout(Exception): pass

  def __init__(self, sec):
    self.sec = sec

  def __enter__(self):
    signal.signal(signal.SIGALRM, self.raise_timeout)
    signal.alarm(self.sec)

  def __exit__(self, *args):
    signal.alarm(0) # disable alarm

  def raise_timeout(self, *args):
    raise Timeout.Timeout()

class OpenKnot(SpaceCurve):
    '''
    Class for holding the vertices of a single line that is assumed to
    be an open curve. This class inherits from
    :class:`~pyknot2.spacecurves.spacecurve.SpaceCurve`, replacing any
    default arguments that assume closed curves, and providing methods
    for statistical analysis of knot invariants on projection and closure.

    All knot invariant methods return the results of a sampling over
    many projections of the knot, unless indicated otherwise.
    '''

    def __init__(self, *args, **kwargs):
        super(OpenKnot, self).__init__(*args, **kwargs)
        self._cached_v2 = {}

    @SpaceCurve.points.setter
    def points(self, points):
        super(OpenKnot, self.__class__).points.fset(self, points)
        self._cached_alexanders = None

    def closing_distance(self):
        '''
        Returns the distance between the first and last points.
        '''
        return n.linalg.norm(self.points[-1] - self.points[0])

    def raw_crossings(self, mode='use_max_jump',
                      virtual_closure=False,
                      recalculate=False, try_cython=False,
                      number_closure_points=1):
        '''
        Calls :meth:`pyknot2.spacecurves.spacecurve.SpaceCurve.raw_crossings`,
        but without including the closing line between the last
        and first points (i.e. setting include_closure=False).
        '''
        if not virtual_closure:
            return super(OpenKnot, self).raw_crossings(mode=mode,
                                                       include_closure=False,
                                                       recalculate=recalculate,
                                                       try_cython=try_cython)
        cs = super(OpenKnot, self).raw_crossings(mode=mode,
                                                 include_closure=True,
                                                 recalculate=recalculate,
                                                 try_cython=try_cython)

        if len(cs) > 0:
            closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-number_closure_points)) |
                                    ((cs[:, 1] > len(self.points)-number_closure_points)))
            indices = closure_cs.flatten()
            for index in indices:
                cs[index, 2] = 0
        return cs

    def __str__(self):
        if self._crossings is not None:
            return '<OpenKnot with {} points, {} crossings>'.format(
                len(self.points), len(self._crossings))
        return '<OpenKnot with {} points>'.format(len(self.points))

    def arclength(self):
        '''Calls :meth:`pyknot2.spacecurves.spacecurve.SpaceCurve.arclength`,
        automatically *not* including the closure.
        '''
        return super(OpenKnot, self).arclength(include_closure=False)

    def smooth(self, repeats=1, window_len=10, window='hanning'):
        '''Calls :meth:`pyknot2.spacecurves.spacecurve.SpaceCurve.smooth`,
        with the periodic argument automatically set to False.
        '''
        super(OpenKnot, self).smooth(repeats=repeats, window_len=window_len,
                                     window=window, periodic=False)

    def _plot_uniform_angles(self, number_of_samples):
        '''
        Plots the projection of the knot at each of the given
        number of samples, approximately evenly distributed on
        the sphere.

        This function is really just for testing.
        '''
        angles = get_rotation_angles(number_of_samples)

        for i, angs in enumerate(angles):
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            fig, ax = k.plot_projection(show=False)
            fig.set_size_inches((2, 2))
            fig.savefig('rotation{:05d}'.format(i))
            fig.close()

    def plot_projections(self, number_of_samples):
        '''
        Plots the projection of the knot at each of the given
        number of samples squared, rotated such that the sample
        direction is vertical.

        The output (and return) is a matplotlib plot with
        number_of_samples x number_of_samples axes.
        '''
        angles = get_rotation_angles(number_of_samples ** 2)

        print('Got angles')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(nrows=number_of_samples,
                                 ncols=number_of_samples)

        print('got axes')

        all_axes = [ax for row in axes for ax in row]

        for i, angs in enumerate(angles):
            self._vprint('i = {} / {}'.format(i, len(angles)))
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            ax = all_axes[i]
            fig, ax = k.plot_projection(fig_ax=(fig, ax), show=False)

        fig.tight_layout()
        fig.show()
        return fig, ax

    def alexander_polynomials(self, number_of_samples=10, radius=None,
                              recalculate=False,
                              zero_centroid=False,
                              optimise_closure=True,
                              random_closure_path=True):
        '''
        Returns a list of Alexander polynomials for the knot, closing
        on a sphere of the given radius, with the given number of sample
        points approximately evenly distributed on the sphere.

        The results are cached by number of samples and radius.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to sample. Defaults to 10.
        optimise_closure: bool
            If True, doesn't really close on a sphere but at infinity.
            This lets the calculation be optimised slightly, and so is the
            default.
        radius : float
            The radius of the sphere on which to close the knot. Defaults
            to None, which picks 10 times the largest Cartesian deviation
            from 0. This is *only* used if optimise_closure=False.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to True.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and alexander
            polynomials.
        '''

        if zero_centroid:
            self.zero_centroid()

        if not recalculate and self._cached_alexanders is not None:
            if (number_of_samples, radius) in self._cached_alexanders:
                return self._cached_alexanders[(number_of_samples,
                                                radius)]
        else:
            self._cached_alexanders = {}
        angles = get_rotation_angles(number_of_samples)

        polys = []

        cache_radius = radius

        if radius is None:
            radius = 100 * n.max(self.points)
            # Not guaranteed to give 10* the real radius, but good enough

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            if random_closure_path:
                mid_point = self.points[0]/2. + self.points[-1]/2.
                mid_point += n.random.random(3)
                k = Knot(n.append(self.points, [mid_point], axis=0), verbose=False)
                number_closure_points = 2
            else:
                k = Knot(self.points, verbose=False)
                number_closure_points = 1
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            if optimise_closure:
                cs = k.raw_crossings()
                if len(cs) > 0:
                    closure_cs = n.argwhere(((cs[:, 0] > len(k.points)-number_closure_points) & (cs[:, 2] < 0.)) |
                                            ((cs[:, 1] > len(k.points)-number_closure_points) & (cs[:, 2] > 0.)))
                    indices = closure_cs.flatten()
                    for index in indices:
                        cs[index, 2:] *= -1
                gc = GaussCode(cs, verbose=self.verbose)
                gc.simplify()
                polys.append([angs[0], angs[1], alexander(gc, simplify=False)])

            else:
                points = k.points
                closure_point = points[-1] + points[0] / 2.
                closure_point[2] = radius
                k.points = n.vstack([points, closure_point])

                polys.append([angs[0], angs[1], k.alexander_polynomial()])

        self._cached_alexanders[
            (number_of_samples, cache_radius)] = n.array(polys)

        return n.array(polys)
    
    def closure_alexander_polynomial(self, theta=0, phi=0, random_closure_path=True):
        '''Returns the Alexander polynomial of the knot, when projected in
        the z plane after rotating the given theta and phi to the
        North pole.

        Parameters
        ----------
        theta : float
            The sphere angle theta
        phi : float
            The sphere angle phi
        '''
        if random_closure_path:
            mid_point = self.points[0]/2. + self.points[-1]/2.
            mid_point += n.random.random(3)
            k = Knot(n.append(self.points, [mid_point], axis=0), verbose=False)
            closure_path_points = 2
        else:
            k = Knot(self.points, verbose=False)
            closure_path_points = 1
        k._apply_matrix(rotate_to_top(theta, phi))

        cs = k.raw_crossings()
        if len(cs) > 0:
            closure_cs = n.argwhere(((cs[:, 0] > len(k.points)-closure_path_points) & (cs[:, 2] < 0.)) |
                                    ((cs[:, 1] > len(k.points)-closure_path_points) & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            for index in indices:
                cs[index, 2:] *= -1
        gc = GaussCode(cs, verbose=False)
        gc.simplify()
        return alexander(gc, simplify=False)

    def alexander_fractions(self, number_of_samples=10, **kwargs):
        '''Returns each of the Alexander polynomials from
        self.alexander_polynomials, with the fraction of that type.
        '''
        polys = self.alexander_polynomials(
            number_of_samples=number_of_samples, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))

        return sorted(fracs, key=lambda j: j[1])

    def _alexander_map_values(self, number_of_samples=10, interpolation=100,
                              **kwargs):
        polys = self.alexander_polynomials(
            number_of_samples=number_of_samples, **kwargs)

        from scipy.interpolate import griddata

        positions = []
        for i, row in enumerate(polys):
            positions.append(gall_peters(row[0], row[1]))
        positions = n.array(positions)

        interpolation_points = n.mgrid[0:2*n.pi:int(1.57*interpolation)*1j,
                                       -2.:2.:interpolation*1j]
        # interpolation_points = n.mgrid[0:2 * n.pi:157j,
        #                        -2.:2.:100j]
        # '''interpolation_points = n.mgrid[-n.pi:n.pi:157j,
        #                        -n.pi/2:n.pi/2:100j]'''

        values = griddata(positions, polys[:, 2],
                          tuple(interpolation_points),
                          method='nearest')

        return positions, values

    def plot_alexander_map(self, number_of_samples=10,
                           scatter_points=False,
                           mode='imshow', interpolation=100,
                           **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different Alexander polynomial in a different colour according
        to a closure on a far away point in this direction.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to close at.
        scatter_points : bool
            If True, plots a dot at each point on the map projection
            where a closure was made.
        mode : str
            'imshow' to plot the pixels of an image, otherwise plots
            filled contours. Defaults to 'imshow'.
        interpolation : int
            The (short) side length of the interpolation grid on which
            the map projection is made. Defaults to 100.
        '''

        positions, values = self._alexander_map_values(
            number_of_samples,
            interpolation=interpolation)

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        '''fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection="mollweide")
        ax.grid(True)'''

        if mode == 'imshow':
            cax = ax.imshow(values.T, cmap='jet', interpolation='none')
            fig.colorbar(cax)
        else:
            ax.contourf(values.T, cmap='jet',
                        levels=[0] + range(3, int(n.max(values) + 1.1), 2))
        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_xlim(-0.5, 1.565*interpolation)
        ax.set_ylim(-0.5, 0.995*interpolation)

        im_positions = positions*25*float(interpolation / 100.)
        im_positions[:, 0] -= 0.5
        im_positions[:, 1] += 49.*(interpolation / 100.) + 0.5
        if scatter_points:
            ax.scatter(im_positions[:, 0], im_positions[:, 1], color='black',
                       alpha=1, s=1)

        fig.tight_layout()
        fig.show()
        return fig, ax

    def plot_alexander_shell(self, number_of_samples=10, radius=None,
                             zero_centroid=False,
                             sphere_radius_factor=2.,
                             opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured by the type of knot obtained by closing on each
        point.

        Parameters are all passed to :meth:`OpenKnot.alexander_polynomials`,
        except opacity and kwargs which are given to mayavi.mesh, and
        sphere_radius_factor which gives the radius of the enclosing
        sphere in terms of the maximum Cartesian distance of any point
        in the line from the origin.
        '''

        self.plot(mode='mayavi')

        positions, values = self._alexander_map_values(
            number_of_samples, radius=None, zero_centroid=False)

        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
        phis = n.linspace(0, 2 * n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor * n.max(self.points)
        zs = r * n.cos(thetas)
        xs = r * n.sin(thetas) * n.cos(phis)
        ys = r * n.sin(thetas) * n.sin(phis)

        import mayavi.mlab as may

        may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)

    def plot_alexander_shell_vispy(self, number_of_samples=10, radius=None,
                                   zero_centroid=False,
                                   sphere_radius_factor=2.,
                                   opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured by the type of knot obtained by closing on each
        point.

        Parameters are all passed to :meth:`OpenKnot.alexander_polynomials`,
        except opacity and kwargs which are given to mayavi.mesh, and
        sphere_radius_factor which gives the radius of the enclosing
        sphere in terms of the maximum Cartesian distance of any point
        in the line from the origin.
        '''

        self.plot(mode='vispy')

        positions, values = self._alexander_map_values(
            number_of_samples, radius=None, zero_centroid=False, **kwargs)

        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi/2.
        phis = n.linspace(0, 2*n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor*n.max(self.points)
        zs = r*n.cos(thetas)
        xs = r*n.sin(thetas)*n.cos(phis)
        ys = r*n.sin(thetas)*n.sin(phis)

        colours = n.zeros((values.shape[0], values.shape[1], 4))
        max_val = n.max(values)
        min_val = n.min(values)
        unique_values = n.unique(colours)
        max_val += (1. + 1./len(unique_values))*(max_val - min_val)
        diff = (max_val - min_val)

        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('jet')
        for i in range(colours.shape[0]):
            for j in range(colours.shape[1]):
                colours[i, j] = cmap(((values[i, j] - min_val) / diff))

        colours[:, :, -1] = opacity

        print('shape is', values.shape)
        from vispy.scene import GridMesh
        from pyknot2.visualise import vispy_canvas
        mesh = GridMesh(xs, ys, zs, colors=colours)
        vispy_canvas.view.add(mesh)

    def virtual_check(self):
        '''
        Takes an open curve and checks (for the default projection) if its
        Gauss code corresponds to a virtual knot or not. Returns a Boolean of
        this information.

        Returns
        -------
        virtual : bool
            True if the Gauss code corresponds to a virtual knot. False
            otherwise.
        '''

        print('Function currently incorrect. Possibly worth correcting')

        # gauss_code = self.gauss_code()._gauss_code[0][:, 0]
        # l = len(gauss_code)
        # total_crossings = l / 2
        # crossing_counter = 1
        # virtual = False
        # seen = set()
        # seen_add = seen.add
        # gauss_code_crossing_numbers = [x for x in gauss_code if not (x in seen or seen_add(x))]
        #
        # for crossing_number in gauss_code_crossing_numbers:
        #     occurences = n.where(gauss_code == crossing_number)[0]
        #     first_occurence = occurences[0]
        #     second_occurence = occurences[1]
        #     crossing_difference = second_occurence - first_occurence
        #
        #     if crossing_difference % 2 == 0:
        #         return True
        # return False

    def virtual_checks(self, number_of_samples=10,
                                        zero_centroid=False):
        '''
        Returns a list of virtual Booleans for the curve with a given number
        if projections taken from directions approximately evenly distributed.
        A value of True corresponds to the projection giving a virtual knot,
        with False returned otherwise.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and virtual Booleans
            (True if virtual, False otherwise)
        '''

        print('Function currently incorrect. Possibly worth correcting')

        # if zero_centroid:
        #     self.zero_centroid()
        #
        # angles = get_rotation_angles(number_of_samples)
        #
        # polys = []
        #
        # print_dist = int(max(1, 3000. / len(self.points)))
        # for i, angs in enumerate(angles):
        #     if i % print_dist == 0:
        #         self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
        #     k = OpenKnot(self.points, verbose=False)
        #     k._apply_matrix(rotate_to_top(*angs))
        #     if zero_centroid:
        #         k.zero_centroid()
        #     isvirtual = k.virtual_check()
        #     polys.append([angs[0], angs[1], isvirtual])
        #
        # return n.array(polys)

    def virtual_fractions(self, number_of_samples=10, **kwargs):
        '''Returns each of the virtual booleans from
        self.virtual.check.projections, with the fraction of each type.
        '''

        print('Function currently incorrect. Possibly worth correcting')

        # polys = self.virtual_checks(
        #     number_of_samples=number_of_samples, **kwargs)
        # alexs = n.round(polys[:, 2]).astype(n.int)
        #
        # fracs = []
        # length = float(len(alexs))
        # for alex in n.unique(alexs):
        #     fracs.append((alex, n.sum(alexs == alex) / length))
        #
        # return sorted(fracs, key=lambda j: j[1])

    def _virtual_map_values(self, number_of_samples=10, **kwargs):

        print('Function currently incorrect. Possibly worth correcting')

        # polys = self.virtual_checks(
        #     number_of_samples=number_of_samples, **kwargs)
        #
        # from scipy.interpolate import griddata
        #
        # positions = []
        # for i, row in enumerate(polys):
        #     positions.append(gall_peters(row[0], row[1]))
        # positions = n.array(positions)
        #
        # interpolation_points = n.mgrid[0:2 * n.pi:157j,
        #                        -2.:2.:100j]
        # values = griddata(positions, polys[:, 2],
        #                   tuple(interpolation_points),
        #                   method='nearest')
        #
        # return positions, values

    def plot_virtual_map(self, number_of_samples=10,
                         scatter_points=False,
                         mode='imshow', **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different virtual Boolean in a different colour according
        to a projection in this direction.
        '''

        print('Function currently incorrect. Possibly worth correcting')

        # positions, values = self._virtual_map_values(number_of_samples)
        #
        # import matplotlib.pyplot as plt
        #
        # fig, ax = plt.subplots()
        #
        # if mode == 'imshow':
        #     cax = ax.imshow(values.T, cmap='jet', interpolation='none')
        #     fig.colorbar(cax)
        # else:
        #     ax.contourf(values.T, cmap='jet',
        #                 levels=[0] + range(3, int(n.max(values) + 1.1), 2))
        # ax.set_xticks([])
        # ax.set_yticks([])
        #
        # ax.set_xlim(-0.5, 156.5)
        # ax.set_ylim(-0.5, 99.5)
        #
        # im_positions = positions * 25
        # im_positions[:, 0] -= 0.5
        # im_positions[:, 1] += 49.5
        # if scatter_points:
        #     ax.scatter(im_positions[:, 0], im_positions[:, 1], color='black',
        #                alpha=1, s=1)
        #
        # fig.tight_layout()
        # fig.show()
        #
        # return fig, ax

    def plot_virtual_shell(self, number_of_samples=10,
                           zero_centroid=False,
                           sphere_radius_factor=2.,
                           opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured according to whether or not the projection from this
        point corresponds to a virtual knot or not.

        Parameters are all passed to
        :meth:`OpenKnot.virtual_checks`, except opacity and kwargs
        which are given to mayavi.mesh, and sphere_radius_factor which gives
        the radius of the enclosing sphere in terms of the maximum Cartesian
        distance of any point in the line from the origin.
        '''

        print('Function currently incorrect. Possibly worth correcting')

        # self.plot()
        #
        # positions, values = self._virtual_map_values(
        #     number_of_samples, zero_centroid=False)
        #
        # thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
        # phis = n.linspace(0, 2 * n.pi, 157)
        #
        # thetas, phis = n.meshgrid(thetas, phis)
        #
        # r = sphere_radius_factor * n.max(self.points)
        # zs = r * n.cos(thetas)
        # xs = r * n.sin(thetas) * n.cos(phis)
        # ys = r * n.sin(thetas) * n.sin(phis)
        #
        # import mayavi.mlab as may
        #
        # may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)

    def self_linking(self, theta=0, phi=0, sphere_closure=False):
        k = OpenKnot(self.points, verbose=False)
        k._apply_matrix(rotate_to_top(theta, phi))

        gausscode = k.gauss_code(virtual_closure=True)
        gausscode.simplify()
        gausscode = gausscode.without_virtual()._gauss_code[0]
        # gausscode = k.gauss_code()._gauss_code[0]
        l = len(gausscode)

        self_linking_counter = 0

        cache = {}

        for index, row in enumerate(gausscode):
            number, over_under, orientation = row
            if number in cache:
                if ((index - cache[number]) % 2) == 0:
                    self_linking_counter += orientation
            else:
                cache[number] = index

        return self_linking_counter

    def self_linkings(self, number_of_samples=10,
                      zero_centroid=False):
        '''
        Returns a list of self linking numbers for the curve with a given
        number of projections taken from directions approximately evenly
        distributed.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and self linking number
        '''
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        polys = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            self_linking = k.self_linking()
            polys.append([angs[0], angs[1], self_linking])

        return n.array(polys)

    def self_linking_fractions(self, number_of_samples=10, **kwargs):
        '''Returns each of the self linking numbers from
        self.virtual.self_link.projections, with the fraction of each type.
        '''
        polys = self.self_linkings(
            number_of_samples=number_of_samples, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]
        
    def _self_linking_map_values(self, number_of_samples=10, **kwargs):
        polys = self.self_linkings(
            number_of_samples=number_of_samples, **kwargs)

        from scipy.interpolate import griddata

        positions = []
        for i, row in enumerate(polys):
            positions.append(gall_peters(row[0], row[1]))
        positions = n.array(positions)

        interpolation_points = n.mgrid[0:2 * n.pi:157j,
                               -2.:2.:100j]
        values = griddata(positions, polys[:, 2],
                          tuple(interpolation_points),
                          method='nearest')

        return positions, values

    def plot_self_linking_map(self, number_of_samples=10,
                           scatter_points=False,
                           mode='imshow', **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different self linking number in a different colour according
        to a projection in this direction.
        '''

        absolute = kwargs.get('absolute', False)
        colour_map = kwargs.get('cmap', 'jet')

        positions, values = self._self_linking_map_values(number_of_samples)

        if absolute:
            values = n.array([abs(value) for value in values])

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        if mode == 'imshow':
            cax = ax.imshow(values.T, cmap=colour_map, interpolation='none')
            fig.colorbar(cax)
        else:
            ax.contourf(values.T, cmap=colour_map,
                        levels=[0] + range(3, int(n.max(values) + 1.1), 2))
        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_xlim(-0.5, 156.5)
        ax.set_ylim(-0.5, 99.5)

        im_positions = positions * 25
        im_positions[:, 0] -= 0.5
        im_positions[:, 1] += 49.5
        if scatter_points:
            ax.scatter(im_positions[:, 0], im_positions[:, 1], color='black',
                       alpha=1, s=1)

        fig.tight_layout()
        fig.show()

        return fig, ax

    def plot_self_linking_shell(self, number_of_samples=10,
                                zero_centroid=False,
                                sphere_radius_factor=2.,
                                opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured by the self linking number obtained by projecting from
        this point.

        Parameters are all passed to 
        :meth:`OpenKnot.virtual_checks`, except opacity and kwargs 
        which are given to mayavi.mesh, and sphere_radius_factor which gives 
        the radius of the enclosing sphere in terms of the maximum Cartesian
        distance of any point in the line from the origin.
        '''

        self.plot()

        positions, values = self._self_linking_map_values(
            number_of_samples, zero_centroid=False)

        absolute = kwargs.get('absolute', False)
        if absolute:
            values = n.array([abs(value) for value in values])


        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
        phis = n.linspace(0, 2 * n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor * n.max(self.points)
        zs = r * n.cos(thetas)
        xs = r * n.sin(thetas) * n.cos(phis)
        ys = r * n.sin(thetas) * n.sin(phis)

        import mayavi.mlab as may

        may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)
        
    def plot_self_linking_shell_vispy(self, number_of_samples=10,
                                      zero_centroid=False,
                                      sphere_radius_factor=2.,
                                      opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured by the self linking number obtained by projecting from
        this point.

        Parameters are all passed to
        :meth:`OpenKnot.virtual_checks`, except opacity and kwargs
        which are given to mayavi.mesh, and sphere_radius_factor which gives
        the radius of the enclosing sphere in terms of the maximum Cartesian
        distance of any point in the line from the origin.
        '''
        self.plot(mode='vispy')

        colour_map = kwargs.get('cmap', 'jet')


        positions, values = self._self_linking_map_values(
            number_of_samples, zero_centroid=False)

        absolute = kwargs.get('absolute', False)
        if absolute:
            values = n.array([abs(value) for value in values])

        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi/2.
        phis = n.linspace(0, 2*n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor*n.max(self.points)
        zs = r*n.cos(thetas)
        xs = r*n.sin(thetas)*n.cos(phis)
        ys = r*n.sin(thetas)*n.sin(phis)

        colours = n.zeros((values.shape[0], values.shape[1], 4))
        max_val = n.max(values)
        min_val = n.min(values)
        unique_values = n.unique(colours)
        #max_val += (1./len(unique_values))*(max_val - min_val)
        diff = (max_val - min_val)

        import matplotlib.pyplot as plt
        cmap = plt.get_cmap(colour_map)
        for i in range(colours.shape[0]):
            for j in range(colours.shape[1]):
                colours[i, j] = cmap(((values[i, j] - min_val) / diff))

        colours[:, :, -1] = opacity

        print('shape is', values.shape)
        from vispy.scene import GridMesh
        from pyknot2.visualise import vispy_canvas
        mesh = GridMesh(xs, ys, zs, colors=colours)
        vispy_canvas.view.add(mesh)

    def average_self_linking(self, number_of_samples=10, **kwargs):
        self_linkings = self.self_linkings(number_of_samples, **kwargs)
        return n.average(self_linkings[:, 2])

    def planar_writhes(self, number_of_samples=10, sphere_closure=False, zero_centroid=False):
        '''
        Returns a list of planar writhes for the curve with a given
        number of projections taken from directions approximately evenly
        distributed.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and self linking number
        '''

        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        writhes = []

        if not sphere_closure:
            print_dist = int(max(1, 3000. / len(self.points)))
            for i, angs in enumerate(angles):
                if i % print_dist == 0:
                    self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
                k = OpenKnot(self.points, verbose=False)
                k._apply_matrix(rotate_to_top(*angs))
                if zero_centroid:
                    k.zero_centroid()
                writhes.append([angs[0], angs[1], k.planar_writhe()])

        else:
            print_dist = int(max(1, 3000. / len(self.points)))
            for i, angs in enumerate(angles):
                if i % print_dist == 0:
                    self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
                open_k = OpenKnot(self.points, verbose=False)
                open_k._apply_matrix(rotate_to_top(*angs))

                mid_point = open_k.points[0]/2. + open_k.points[-1]/2.
                mid_point += n.random.random(3)

                k = Knot(open_k.points, verbose=False)
                number_closure_points = 2
                cs = k.raw_crossings()
                if len(cs) > 0:
                    # close with over crossings
                    closure_cs = n.argwhere(((cs[:, 1] > len(open_k.points)-number_closure_points)
                                            & (cs[:, 2] < 0.)) |
                                            ((cs[:, 0] > len(open_k.points)-number_closure_points)
                                            & (cs[:, 2] > 0.)))
                    indices = closure_cs.flatten()
                    over_cs = copy.deepcopy(cs)
                    for index in indices:
                        over_cs[index, 2:] *= -1
                    closed_over_gc = GaussCode(over_cs, verbose=open_k.verbose)
                    writhes.append([angs[0], angs[1], n.sum(closed_over_gc._gauss_code[0][:, 2])/2.])
                else:
                    writhes.append([angs[0], angs[1], 0])

        return n.array(writhes)

    def average_planar_writhe(self, number_of_samples=10, sphere_closure=False, **kwargs):
        p_writhes = self.planar_writhes(number_of_samples, sphere_closure, **kwargs)
        return n.average(p_writhes[:, 2])

    def alexander_polynomials_multiroots(self, number_of_samples=10,
                                         radius=None,
                                         zero_centroid=False):
        '''
        Returns a list of Alexander polynomials for the knot, closing
        on a sphere of the given radius, with the given number of sample
        points approximately evenly distributed on the sphere. The
        Alexander polynomials are found at three different roots (2, 3
        and 4) and the knot types corresponding to these roots are
        returned also.

        The results are cached by number of samples and radius.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to sample. Defaults to 10.
        radius : float
            The radius of the sphere on which to close the knot. Defaults
            to None, which picks 10 times the largest Cartesian deviation
            from 0.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to True.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and alexander
            polynomials.
        '''
        from pyknot2.catalogue.identify import from_invariants
        from pyknot2.catalogue.database import Knot as dbknot
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        polys = []

        if radius is None:
            radius = 100 * n.max(self.points)
            # Not guaranteed to give 10* the real radius, but good enough

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = Knot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            points = k.points
            closure_point = points[-1] + points[0] / 2.
            closure_point[2] = radius
            k.points = n.vstack([points, closure_point])
            root_at_two = k.alexander_at_root(2)
            root_at_three = k.alexander_at_root(3)
            root_at_four = k.alexander_at_root(4)
            k_gauss_code = k.gauss_code()
            k_gauss_code.simplify()
            max_crossings = len(k_gauss_code)
            if max_crossings > 17:
                max_crossings = 18
            knot_type = from_invariants(determinant=root_at_two, alex_imag_3=root_at_three,
                                        alex_imag_4=root_at_four,
                                        other=[dbknot.min_crossings <= max_crossings])
            knot_type = knot_db_to_string(knot_type)
            polys.append([angs[0], angs[1], root_at_two, root_at_three, root_at_four, max_crossings,
                          knot_type])

        return polys

    def alexander_multiroots_fractions(self, number_of_samples=10, **kwargs):

        '''
        Returns each of the knot types from
        self.alexander_polynomials_multiroots, with the fraction of that type.
        '''

        knot_info = self.alexander_polynomials_multiroots(
            number_of_samples, **kwargs)
        knot_types = []
        for closures in knot_info:
            knot_types.append(closures[6])
        knot_types = [item for sublist in knot_types for item in sublist] #flattens list
        knot_frequency = Counter(knot_types)
        common = knot_frequency.most_common()
        list_common = [list(elem) for elem in common]
        list_common_fractions = [[elem[0],elem[1]/float(number_of_samples)] for elem in list_common]
        return list_common_fractions

    def generalised_alexander_symbolic(self):
        '''
        Returns the generalised Alexander polynomial for the default projection
        of the open knot
        '''

        import sympy as sym

        gauss_code = self.gauss_code()
        gauss_code_crossings = gauss_code._gauss_code[0][:, 0]
        gauss_code_over_under = gauss_code._gauss_code[0][:, 1]
        gauss_code_orientations = gauss_code._gauss_code[0][:, 2]
        seen = set()
        seen_add = seen.add
        gauss_code_crossing_numbers = [x for x in gauss_code_crossings if not (x in seen or seen_add(x))]
        x = sym.var('x')
        y = sym.var('y')
        m_plus = sym.Matrix([[1 - x, -y], [-x * y**-1, 0]])
        m_minus = sym.Matrix([[0, -x**-1 * y], [-y**-1, 1 - x**-1]])
        num_crossings = len(self.gauss_code())
        matrix = sym.zeros(2 * num_crossings, 2 * num_crossings)
        permutation_matrix = sym.zeros(2 * num_crossings, 2 * num_crossings)

        arc_labels = [0]*len(gauss_code_crossings)
        for i in range(len(gauss_code_crossings)):
            arc_labels[i] = [0] * 4

        counter = 0
        for crossing_number in gauss_code_crossing_numbers:
            occurrences = n.where(gauss_code_crossings == crossing_number)[0]
            if gauss_code_orientations[occurrences[0]] == 1:
                m = m_plus
            else:
                m = m_minus
            for i in [0, 1]:
                for j in [0, 1]:
                    matrix[counter*2 + i, counter*2 + j] = m[i, j]
            counter += 1

        for i in range(len(gauss_code_crossings)):
            arc_labels[i][0] = gauss_code_crossings[i]
            arc_labels[i][1] = (gauss_code_orientations[i] *
                                gauss_code_over_under[i])  # -1 = r, +1 = l
            arc_labels[i-1][2] = gauss_code_crossings[i]
            arc_labels[i-1][3] = (-1 * gauss_code_orientations[i] *
                                  gauss_code_over_under[i])  # -1 = r, +1 = l

        counter = 1
        for crossing_number in gauss_code_crossing_numbers:
            for i in range(len(gauss_code_crossings)):
                if arc_labels[i][0] == crossing_number:
                    arc_labels[i][0] = counter
                if arc_labels[i][2] == crossing_number:
                    arc_labels[i][2] = counter
            counter += 1

        for i in range(len(arc_labels)):
            if arc_labels[i][1] < 0:
                if arc_labels[i][3] < 0:
                    # bottom right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-1] = 1
                else:
                    # bottom left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-1] = 1
            else:
                if arc_labels[i][3] < 0:
                    # upper right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-2] = 1
                else:
                    # upper left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-2] = 1

        writhe = sum(gauss_code_orientations)/2

        return (-1)**writhe * ((matrix - permutation_matrix).det())

    def generalised_alexander(self, root_s=3, root_t=4):
        '''
        Returns the generalised Alexander polynomial for the default projection
        at the given roots of unity, i.e. evaluated at exp(2 pi I / root). By
        default root_s = 3 and root_t = 4

        Parameters
        ----------
        root_s and root_t: int
            The roots of unity to use, i.e. evaluating at exp(2 pi I / root).
            If this is iterable, this method returns a list of the results
            at every value of that iterable.
        '''

        root_unity_s = n.exp(2 * n.pi * 1.j / root_s)
        root_unity_s = root_unity_s.real.round(3) + root_unity_s.imag.round(3)*1.j
        root_unity_t = n.exp(2 * n.pi * 1.j / root_t)
        root_unity_t = root_unity_t.real.round(3) + root_unity_t.imag.round(3)*1.j

        gauss_code = self.gauss_code()
        gauss_code_crossings = gauss_code._gauss_code[0][:, 0]
        gauss_code_over_under = gauss_code._gauss_code[0][:, 1]
        gauss_code_orientations = gauss_code._gauss_code[0][:, 2]
        seen = set()
        seen_add = seen.add
        gauss_code_crossing_numbers = [x for x in gauss_code_crossings if not (x in seen or seen_add(x))]
        s = root_unity_s
        t = root_unity_t
        x = s*t
        y = -t

        m_plus = n.matrix([[1 - x, -y], [-x * y**-1, 0]])
        m_minus = n.matrix([[0, -x**-1 * y], [-y**-1, 1 - x**-1]])
        num_crossings = len(self.gauss_code())
        matrix = n.zeros((2 * num_crossings, 2 * num_crossings)) + 0.j
        permutation_matrix = n.zeros((2 * num_crossings, 2 * num_crossings)) + 0.j

        arc_labels = [0]*len(gauss_code_crossings)
        for i in range(len(gauss_code_crossings)):
            arc_labels[i] = [0] * 4

        counter = 0
        for crossing_number in gauss_code_crossing_numbers:
            occurrences = n.where(gauss_code_crossings == crossing_number)[0]
            if gauss_code_orientations[occurrences[0]] == 1:
                m = m_plus
            else:
                m = m_minus
            for i in [0, 1]:
                for j in [0, 1]:
                    matrix[counter*2 + i, counter*2 + j] = m[i, j]
            counter += 1

        for i in range(len(gauss_code_crossings)):
            arc_labels[i][0] = gauss_code_crossings[i]
            arc_labels[i][1] = (gauss_code_orientations[i] *
                                gauss_code_over_under[i])  # -1 = r, +1 = l
            arc_labels[i-1][2] = gauss_code_crossings[i]
            arc_labels[i-1][3] = (-1 * gauss_code_orientations[i] *
                                  gauss_code_over_under[i])  # -1 = r, +1 = l

        counter = 1
        for crossing_number in gauss_code_crossing_numbers:
            for i in range(len(gauss_code_crossings)):
                if arc_labels[i][0] == crossing_number:
                    arc_labels[i][0] = counter
                if arc_labels[i][2] == crossing_number:
                    arc_labels[i][2] = counter
            counter += 1

        for i in range(len(arc_labels)):
            if arc_labels[i][1] < 0:
                if arc_labels[i][3] < 0:
                    # bottom right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-1] = 1
                else:
                    # bottom left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-1] = 1
            else:
                if arc_labels[i][3] < 0:
                    # upper right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-2] = 1
                else:
                    # upper left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-2] = 1

        writhe = sum(gauss_code_orientations)/2

        return n.linalg.det(int((-1)**writhe) * (matrix - permutation_matrix))

    def generalised_alexander_at_root(self, root_s, root_t, round=True):
        '''
        Returns the generalised Alexander polynomial for the default projection
        at the given roots of unity, i.e. evaluated at exp(2 pi I / root).

        The result returned is the absolute value.

        Combinations of roots that work particularly well include:
        [root_s, root_t] = [2, 3], [2, 4] and [3,4]

        Parameters
        ----------
        root_s and root_t : int
            The root of unity to use, i.e. evaluating at exp(2 pi I / root).
            If this is iterable, this method returns a list of the results
            at every value of that iterable.
        round : bool
            If True, the result will be rounded
            to the nearest integer for convenience, and returned as an
            integer type.
        '''

        value = self.generalised_alexander(root_s, root_t)

        if round:
            return int(n.round(n.abs(value)))
        else:
            return n.abs(value)

    def generalised_alexanders(self, number_of_samples=10,
                      root_s=3, root_t=4, zero_centroid=False):
        '''
        Returns a list of generalised alexander polynomials for the curve
        with a given number of projections taken from directions
        approximately evenly distributed.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        root_s and root_t : int
            The root of unity to use, i.e. evaluating at exp(2 pi I / root).
            If this is iterable, this method returns a list of the results
            at every value of that iterable.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and generalised alexander
            polynomial
        '''
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        polys = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            generalised_alexander = k.generalised_alexander_at_root(root_s, root_t)
            polys.append([angs[0], angs[1], generalised_alexander])

        return n.array(polys)

    def generalised_alexanders_fractions(
            self, number_of_samples=10, root_s=3, root_t=4, **kwargs):
        '''Returns each of the self linking numbers from
        self.virtual.self_link.projections, with the fraction of each type.
        '''
        polys = self.generalised_alexanders(
                    number_of_samples, root_s, root_t, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def generalised_alexander_multiroot(self, round = True):
        '''
        Returns the generalised alexander polynomial of the knot for multiple
        combinations of roots of unity. Specifically, (2,3), (2, 4) and (3, 4)
        '''
        first_root = self.generalised_alexander_at_root(2, 3)
        second_root = self.generalised_alexander_at_root(2, 4)
        third_root = self.generalised_alexander_at_root(3, 4)

        return [first_root, second_root, third_root]

    def generalised_alexander_multiroots(self, number_of_samples=10,
                      zero_centroid=False):
        '''
        Returns a list of generalised alexander polynomials for the curve
        with a given number of projections taken from directions
        approximately evenly distributed.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        root_s and root_t : int
            The root of unity to use, i.e. evaluating at exp(2 pi I / root).
            If this is iterable, this method returns a list of the results
            at every value of that iterable.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 3 array of angles and generalised alexander
            polynomial
        '''
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        polys = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            multiroot = k.generalised_alexander_multiroot()
            info = [angs[0], angs[1], multiroot[0], multiroot[1], multiroot[2]]
            polys.append(info)

        return n.array(polys)

    def generalised_alexander_multiroots_fractions(
            self, number_of_samples=10, **kwargs):
        '''Returns each of the generalised alexander for multiple roots from
        self.virtual.generalise_alexander_multiroots, with the fraction of
        each type.
        '''
        polys = self.generalised_alexander_multiroots(
                    number_of_samples, **kwargs)
        alexs = polys[:, 2:5].astype(n.int)

        unique_alexs = [alexs[0]]
        for rows in alexs:
            unique_check = True
            for entries in unique_alexs:
                if (entries == rows).all():
                    unique_check = False
            if unique_check:
                unique_alexs.append(rows)

        fracs = []
        length = float(len(alexs))

        for uniques in unique_alexs:
            occurences = 0
            for alex in alexs:
                if (alex == uniques).all():
                    occurences += 1
            fracs.append((list(uniques), occurences / length))

        # for alex in unique_alexs:
        #    fracs.append((alex, n.sum(alexs == alex, axis=0)[0] / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def jones_polynomial(self, root=-1, simplify=True, random_closure_path=True):
        '''
        Calculates the jones polynomial for a single projection of the open curve
        using the given root
        '''
        import pyknot2.representations.planardiagram as pdiag
        if simplify:
            if random_closure_path:
                mid_point = self.points[0]/2. + self.points[-1]/2.
                mid_point += n.random.random(3)
                k = OpenKnot(n.append(self.points, [mid_point], axis=0), verbose=False)
                gc = k.gauss_code(recalculate=True, virtual_closure=True, number_closure_points=2)
            else:
                gc = self.gauss_code(recalculate=True, virtual_closure=True)
            gc.simplify()
            gc = gc.without_virtual()
        else:
            gc = self.gauss_code()
        diagram = pdiag.PlanarDiagram(gc)
        return diagram.jones_optimised(root_of_unity=root)

    def jones_polynomials(self, number_of_samples=10, root=-1, simplify=True, **kwargs):
        '''
        calculates the jones polynomial over multiple projections
        '''
        angles = get_rotation_angles(number_of_samples)

        polys = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            jones = k.jones_polynomial(root, simplify, **kwargs)
            polys.append([angs[0], angs[1], jones])

        return n.array(polys)

    def jones_polynomials_fractions(self, number_of_samples=10, root=-1, simplify=True, **kwargs):
        '''
        Returns each of the Jones polynomials with the fraction of each.
        '''
        polys = self.jones_polynomials(
            number_of_samples=number_of_samples, root=root, simplify=simplify, **kwargs)
        joneses = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(joneses))
        for jones in n.unique(joneses):
            fracs.append((jones, n.sum(joneses == jones) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def projection_invariant(self, theta=0, phi=0, simplify=True, legacy=True, vassiliev=False,
                             random_closure_path=True):
        '''
        First calculates the generalised Alexander polynomial of the curve. If 0,
        then calculates the classical Alexander determinant and vassiliev invariant
        of degree 2.

        Returned is a list of invariants in the following order:
        [gen_alex(2,3), gen_alex(2,4), gen_alex(3,4), alex(-1), vassiliev_2]
        Returns 0 for classical invariants if projection is virtual
        '''

        open_k = OpenKnot(self.points, verbose=False)
        open_k._apply_matrix(rotate_to_top(theta, phi))

        # if simplify:
        #     gc = open_k.gauss_code(recalculate=True, virtual_closure=True)
        #     gc.simplify()
        #     gc = gc.without_virtual()
        # else:
        #     gc = open_k.gauss_code()

        if not simplify and random_closure_path:
            print('Always simplify with random closure path.\nSimplify set to True')
            simplify = True

        if simplify:
            if random_closure_path:
                mid_point = open_k.points[0]/2. + open_k.points[-1]/2.
                mid_point += n.random.random(3)
                open_k = OpenKnot(n.append(open_k.points, [mid_point], axis=0), verbose=False)
                gc = open_k.gauss_code(recalculate=True, virtual_closure=True, number_closure_points=2)
            else:
                gc = open_k.gauss_code(recalculate=True, virtual_closure=True)
            gc.simplify()
            gc = gc.without_virtual()
        else:
            gc = open_k.gauss_code()

        if len(gc) < 2:
            if not legacy:
                if vassiliev:
                    return [0, 0, 0, 1, 0]
                else:
                    return [0, 0, 0, 1, 1, 1]
            else:
                if vassiliev:
                    return [0, 0, 0, 1, 0, 1, 0]
                else:
                    return [0, 0, 0, 1, 1, 1, 1, 1, 1]

        #if legacy:
        k = Knot(open_k.points, verbose=False)
        if random_closure_path:
            number_closure_points = 2
        else:
            number_closure_points = 1
        cs = k.raw_crossings()
        if len(cs) > 0:
            # close with under crossings
            closure_cs = n.argwhere(((cs[:, 0] > len(open_k.points)-number_closure_points)
                                    & (cs[:, 2] < 0.)) |
                                    ((cs[:, 1] > len(open_k.points)-number_closure_points)
                                    & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            under_cs = copy.deepcopy(cs)
            for index in indices:
                under_cs[index, 2:] *= -1

        # close with over crossings
            closure_cs = n.argwhere(((cs[:, 1] > len(open_k.points)-number_closure_points)
                                    & (cs[:, 2] < 0.)) |
                                    ((cs[:, 0] > len(open_k.points)-number_closure_points)
                                    & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            over_cs = copy.deepcopy(cs)
            for index in indices:
                over_cs[index, 2:] *= -1
        else:
            if legacy:
                if vassiliev:
                    return [0, 0, 0, 1, 0, 1, 0]
                else:
                    return [0, 0, 0, 1, 1, 1, 1, 1, 1]
            else:
                if vassiliev:
                    return [0, 0, 0, 1, 0]
                else:
                    return [0, 0, 0, 1, 1, 1]
        closed_over_gc = GaussCode(over_cs, verbose=open_k.verbose)
        closed_under_gc = GaussCode(under_cs, verbose=open_k.verbose)

        # # Remove closing crossings to calculate self linking
        # xs, ys = n.where(cs[:, :2] > len(k.points) - 1)
        # keeps = n.ones(len(cs), dtype=n.bool)
        # for x in xs:
        #     keeps[x] = False
        # cs = cs[keeps]
        # open_gc = GaussCode(cs)

        gen_alex_1 = generalised_alexander(gc, root_s=2, root_t=3)
        gen_alex_2 = generalised_alexander(gc, root_s=2, root_t=4)
        gen_alex_3 = generalised_alexander(gc, root_s=3, root_t=4)
        if [gen_alex_1, gen_alex_2, gen_alex_3] != [0, 0, 0]: # or open_k.self_linking():  # self linking doesn't seem to detect anything more than gen_alex, consider removing
            if not legacy:
                if vassiliev:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, 0, 0]
                else:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, 0, 0, 0]
            else:
                if vassiliev:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, 0, 0, 0, 0]
                else:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, 0, 0, 0, 0, 0, 0]
        else:
            #alex_closed = alexander(closed_gc, -1)
            alex_over = alexander(closed_over_gc, -1)
            #alex_closed = int(n.round(n.abs(alex_closed)))
            alex_over = int(n.round(n.abs(alex_over)))

            if vassiliev:
                #v2_closed = vassiliev_degree_2(closed_gc)
                v2_over = vassiliev_degree_2(closed_over_gc)
                #v2_closed = int(n.round(n.abs(v2_closed)))
                v2_over = int(n.round(n.abs(v2_over)))
                #return [gen_alex_1, gen_alex_2, gen_alex_3, alex_closed, v2_closed, alex_open, v2_open]
                if not legacy:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, alex_over, v2_over]
                else:
                    alex_under = alexander(closed_under_gc, -1)
                    alex_under = int(n.round(n.abs(alex_under)))
                    v2_under = vassiliev_degree_2(closed_under_gc)
                    v2_under = int(n.round(n.abs(v2_under)))
                    return [gen_alex_1, gen_alex_2, gen_alex_3, alex_over, v2_over, alex_under, v2_under]

            else:
                root_2 = n.exp(2 * n.pi * 1.j / 3.)
                root_2 = root_2.real.round(3) + root_2.imag.round(3)*1.j
                root_3 = n.exp(2 * n.pi * 1.j / 4.)
                root_3 = root_3.real.round(3) + root_3.imag.round(3)*1.j

                alex_over_2 = alexander(gc, root_2)
                alex_over_2 = int(n.round(n.abs(alex_over_2)))
                alex_over_3 = alexander(gc, root_3)
                alex_over_3 = int(n.round(n.abs(alex_over_3)))

                if not legacy:
                    return [gen_alex_1, gen_alex_2, gen_alex_3, alex_over, alex_over_2, alex_over_3]
                else:
                    alex_under = alexander(closed_under_gc, -1)
                    alex_under = int(n.round(n.abs(alex_under)))
                    alex_under_2 = alexander(closed_under_gc, root_2)
                    alex_under_2 = int(n.round(n.abs(alex_under_2)))
                    alex_under_3 = alexander(closed_under_gc, root_3)
                    alex_under_3 = int(n.round(n.abs(alex_under_3)))
                    return [gen_alex_1, gen_alex_2, gen_alex_3,
                            alex_over, alex_over_2, alex_over_3,
                            alex_under, alex_under_2, alex_under_3]

    def projection_invariants(self, number_of_samples=10,
                              zero_centroid=False, legacy=True, vassiliev=False, **kwargs):
        '''
        Returns a list of projection based invariants for the curve
        with a given number of projections taken from directions
        approximately evenly distributed.

        Parameters
        ----------
        number_of_samples : int
            The number of points on the sphere to project from. Defaults to 10.
        zero_centroid : bool
            Whether to first move the average position of vertices to
            (0, 0, 0). Defaults to False.

        Returns
        -------
        : ndarray
            A number_of_samples by 7 array of angles and projection invariants.
            invariants appear in order: gen_alex(2,3), gen_alex(2,4), gen_alex(3,4),
            alex(-1), vassiliev_2
        '''
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        full_invariants = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            invariants = k.projection_invariant(legacy=legacy, vassiliev=vassiliev, **kwargs)
            if not legacy:
                if vassiliev:
                    full_invariants.append([angs[0], angs[1], invariants[0], invariants[1],
                                           invariants[2], invariants[3], invariants[4]])
                else:
                    full_invariants.append([angs[0], angs[1], invariants[0], invariants[1],
                                           invariants[2], invariants[3], invariants[4], invariants[5]])
            else:
                if vassiliev:
                    full_invariants.append([angs[0], angs[1], invariants[0], invariants[1],
                                           invariants[2], invariants[3], invariants[4], invariants[5],
                                           invariants[6]])
                else:
                    full_invariants.append([angs[0], angs[1], invariants[0], invariants[1],
                                           invariants[2], invariants[3], invariants[4], invariants[5],
                                           invariants[6], invariants[7], invariants[8]])

        return n.array(full_invariants)

    def projection_invariants_fractions(self, number_of_samples=10, summary=False, legacy=True, vassiliev=False,
                                        **kwargs):
        '''
        Returns each of the projection invariants from
        self.virtual.projections_invariants, with the fraction of each type. If
        summary=True also returns [unknot %, classical %, virtual %],
        [most common classical invariant, %] and
        [most common virtual invariant, %]
        '''
        if summary and not legacy:
            legacy = True
            print('legacy automatically set to True as summary == True')
        polys = self.projection_invariants(number_of_samples, legacy=legacy, vassiliev=vassiliev, **kwargs)
        full_invariants = polys[:, 2:].astype(n.int)
        # if not legacy:
        #     closure_invariants = polys[:, 2:5].astype(n.int)
        # else:
        #     closure_invariants = polys[:, 2:7].astype(n.int)

        unique_invariants = [full_invariants[0]]
        for rows in full_invariants:
            unique_check = True
            for entries in unique_invariants:
                if (entries == rows).all():
                    unique_check = False
            if unique_check:
                unique_invariants.append(rows)

        fracs = []
        length = float(len(full_invariants))

        for uniques in unique_invariants:
            occurences = 0
            for alex in full_invariants:
                if (alex == uniques).all():
                    occurences += 1
            fracs.append((list(uniques), occurences / length))

        full_invariant_fracs = sorted(fracs, key=lambda j: j[1])

        if summary:
        #     unique_invariants = [closure_invariants[0]]
        #     for rows in closure_invariants:
        #         unique_check = True
        #         for entries in unique_invariants:
        #             if (entries == rows).all():
        #                 unique_check = False
        #         if unique_check:
        #             unique_invariants.append(rows)
        #
        #     fracs = []
        #     length = float(len(full_invariants))
        #
        #     for uniques in unique_invariants:
        #         occurences = 0
        #         for invariant in full_invariants:
        #             if (invariant == uniques).all():
        #                 occurences += 1
        #         fracs.append((list(uniques), occurences / length))
        #
        #     closure_invariant_fracs = sorted(fracs, key=lambda j: j[1])

            classical_number = 0.0
            virtual_number = 0.0
            unknot_number = 0.0
            highest_classical = 0
            if vassiliev:
                classical_mode = [None, None]
            else:
                classical_mode = [None, None, None]
            highest_virtual = 0
            virtual_mode = [None, None, None]

            for projection_info in full_invariant_fracs:
                if (projection_info[0][3:] == [1, 0, 1, 0] and vassiliev) or \
                   projection_info[0][3:] == [1, 1, 1, 1, 1, 1]:
                    unknot_number += projection_info[1]
                elif (vassiliev and (projection_info[0][3:] == [0, 0, 0, 0] or
                      projection_info[0][3:5] != projection_info[0][5:7])) or \
                        (not vassiliev and (projection_info[0][3:] == [0, 0, 0, 0, 0, 0] or
                         projection_info[0][3:6] != projection_info[0][6:9])):  # if this is the case, the projection is virtual
                    virtual_number += projection_info[1]
                    if projection_info[1] >= highest_virtual:
                        virtual_mode = projection_info[0]
                        highest_virtual = projection_info[1]
                else:
                    classical_number += projection_info[1]
                    if projection_info[1] >= highest_classical:
                        classical_mode = projection_info[0]
                        highest_classical = projection_info[1]
                #
                # elif projection_info[0][0:3] == [0, 0, 0]:
                #     classical_number += projection_info[1]
                #     if projection_info[1] >= highest_classical:
                #         classical_mode = projection_info[0]
                #         highest_classical = projection_info[1]
                # else:
                #     virtual_number += projection_info[1]
                #     if projection_info[1] >= highest_virtual:
                #         virtual_mode = projection_info[0]
                #         highest_virtual = projection_info[1]

            return (full_invariant_fracs, [unknot_number, classical_number, virtual_number],
                    [classical_mode, highest_classical], [virtual_mode, highest_virtual])
        else:
            return full_invariant_fracs

    def jones_for_distinguishing(self, root=-1, simplify=True, theta=0, phi=0, legacy=True, timeout=True,
                                 vassiliev=False):
        '''
        First runs generalised_alexander_multiroot, and if the result is
        [3, 4, 5], runs Jones. Otherwise returns 0.

        The reason for this is that virtual knots 2.1, 3.2 and 4.94 all
        have the same generalised Alexander polynomial, but are distinguished
        by Jones. Jones at -1 returns 2, 4 and 5 respectively. It is possible
        that more complicated virtual knots have a Jones at -1 of 0. If this
        becomes an issue, this function will have to be modified (already
        modified from -1 at 1)
        '''

        k = OpenKnot(self.points, verbose=False)
        k._apply_matrix(rotate_to_top(theta, phi))

        alexander_check = k.projection_invariant(legacy=legacy, vassiliev=vassiliev)
        if timeout:
            try:
                with Timeout(300):
                    if not legacy:
                        if alexander_check[:3] == [3, 4, 5]:
                            if vassiliev:
                                return [k.jones_polynomial(root, simplify), '-', '-', '-', '-']
                                #return [k.jones_polynomial(root, simplify), 0, 0, 0, 0]
                            else:
                                return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-']
                        else:
                            return alexander_check
                    else:
                        if alexander_check[:3] == [3, 4, 5]:
                            if vassiliev:
                                return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-', '-']
                                #return [k.jones_polynomial(root, simplify), 0, 0, 0, 0, 0, 0]
                            else:
                                return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-', '-', '-', '-']
                        else:
                            return alexander_check
            except Timeout.Timeout:
                return alexander_check
        else:
            if not legacy:
                if alexander_check[:3] == [3, 4, 5]:
                    if vassiliev:
                        return [k.jones_polynomial(root, simplify), '-', '-', '-', '-']
                        #return [k.jones_polynomial(root, simplify), 0, 0, 0, 0]
                    else:
                        return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-']
                else:
                    return alexander_check
            else:
                if alexander_check[:3] == [3, 4, 5]:
                    if vassiliev:
                        return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-', '-']
                        #return [k.jones_polynomial(root, simplify), 0, 0, 0, 0, 0, 0]
                    else:
                        return [k.jones_polynomial(root, simplify), '-', '-', '-', '-', '-', '-', '-', '-']
                else:
                    return alexander_check

    def jones_for_distinguishings(self, number_of_samples=10, root=-1, simplify=True, legacy=True, vassiliev=False,
                                  **kwargs):
        '''
        calculates the jones_for_distinguishing over multiple projections
        '''
        angles = get_rotation_angles(number_of_samples)

        polys = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            jones = k.jones_for_distinguishing(root, simplify, legacy=legacy, timeout=True, vassiliev=vassiliev)
            polys.append(n.concatenate((angs, jones)))
        return n.array(polys)

    def jones_for_distinguishings_fractions(self, number_of_samples=10, root=-1, simplify=True, legacy=True,
                                            vassiliev=False, summary=False, **kwargs):
        '''
        Returns each of the jones_for_distinguishing with the fraction of each.
        '''

        polys = self.jones_for_distinguishings(
            number_of_samples=number_of_samples, root=root, simplify=simplify, legacy=legacy, vassiliev=vassiliev,
            **kwargs)

        # this is for formatting. makes sure every entry is given as strings instead of floats and the entire array
        # of polys exists as one array, rather than an array of separate, dtype('S32') arrays
        if len(polys) > 1:
            polys = [n.concatenate((poly[:2], [str(number) for number in poly[2:]])) for poly in polys]
            new_polys = [polys[0]]
            for poly in polys[1:]:
                new_polys = n.concatenate((new_polys, [poly]))
            polys = new_polys
        else:
            polys = n.array([[str(x) for x in polys[0]]])

        #if not legacy:
        #    joneses = polys[:, 2:7]
        #else:
        #    joneses = polys[:, 2:9]
        joneses = polys[:, 2:]

        unique_joneses = [joneses[0]]
        for rows in joneses:
            unique_check = True
            for entries in unique_joneses:
                if (entries == rows).all():
                    unique_check = False
            if unique_check:
                unique_joneses.append(rows)

        fracs = []
        length = float(len(joneses))

        for uniques in unique_joneses:
            occurences = 0
            for jones in joneses:
                if (jones == uniques).all():
                    occurences += 1
            fracs.append((list(uniques), occurences / length))

        # for alex in unique_alexs:
        #    fracs.append((alex, n.sum(alexs == alex, axis=0)[0] / length))
        # fracs = n.array(fracs)

        full_invariant_fracs = sorted(fracs, key=lambda j: j[1])

        if summary:
            classical_number = 0.0
            virtual_number = 0.0
            unknot_number = 0.0
            highest_classical = 0
            if vassiliev:
                classical_mode = [None, None]
            else:
                classical_mode = [None, None, None]
            highest_virtual = 0
            virtual_mode = [None, None, None]

            for projection_info in full_invariant_fracs:
                if (projection_info[0][3:] == ['1.0', '0.0', '1.0', '0.0'] and vassiliev) or \
                   projection_info[0][3:] == ['1.0', '1.0', '1.0', '1.0', '1.0', '1.0']:
                    unknot_number += projection_info[1]
                elif (vassiliev and (projection_info[0][3:] == ['0.0', '0.0', '0.0', '0.0'] or
                      projection_info[0][3:5] != projection_info[0][5:7])) or \
                        (not vassiliev and (projection_info[0][3:] == ['0.0', '0.0', '0.0', '0.0', '0.0', '0.0'] or
                         projection_info[0][3:6] != projection_info[0][6:9])) or \
                        (projection_info[0][3] == '-'):  # if this is the case, the projection is virtual
                    virtual_number += projection_info[1]
                    if projection_info[1] >= highest_virtual:
                        virtual_mode = projection_info[0]
                        highest_virtual = projection_info[1]
                else:
                    classical_number += projection_info[1]
                    if projection_info[1] >= highest_classical:
                        classical_mode = projection_info[0]
                        highest_classical = projection_info[1]
                #
                # elif projection_info[0][0:3] == [0, 0, 0]:
                #     classical_number += projection_info[1]
                #     if projection_info[1] >= highest_classical:
                #         classical_mode = projection_info[0]
                #         highest_classical = projection_info[1]
                # else:
                #     virtual_number += projection_info[1]
                #     if projection_info[1] >= highest_virtual:
                #         virtual_mode = projection_info[0]
                #         highest_virtual = projection_info[1]

            return (full_invariant_fracs, [unknot_number, classical_number, virtual_number],
                    [classical_mode, highest_classical], [virtual_mode, highest_virtual])
        else:
            return full_invariant_fracs
        #return fracs[n.argsort(fracs[:, 1])]

    def closure_invariant(self, theta=0, phi=0, vassiliev=False, random_closure_path=True):
        '''
        '''
        if random_closure_path:
            mid_point = self.points[0]/2. + self.points[-1]/2.
            mid_point += n.random.random(3)
            k = Knot(n.append(self.points, [mid_point], axis=0), verbose=False)
            closure_path_points = 2
        else:
            k = Knot(self.points, verbose=False)
            closure_path_points = 1
        k._apply_matrix(rotate_to_top(theta, phi))

        cs = k.raw_crossings()
        if len(cs) > 0:
            closure_cs = n.argwhere(((cs[:, 0] > len(k.points)-closure_path_points) & (cs[:, 2] < 0.)) |
                                    ((cs[:, 1] > len(k.points)-closure_path_points) & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            for index in indices:
                cs[index, 2:] *= -1
        closed_gc = GaussCode(cs, verbose=False)

        # k = Knot(self.points, verbose=False)
        # k._apply_matrix(rotate_to_top(theta, phi))
        #
        # cs = k.raw_crossings()
        # if len(cs) > 0:
        #     closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-1)
        #                             & (cs[:, 2] < 0.)) |
        #                             ((cs[:, 1] > len(self.points)-1)
        #                             & (cs[:, 2] > 0.)))
        #     indices = closure_cs.flatten()
        #     for index in indices:
        #         cs[index, 2:] *= -1
        # closed_gc = GaussCode(cs, verbose=self.verbose)

        try:
            closed_gc.simplify()
        except IndexError:
            print('Could not simplify closed Gauss code\n')

        if len(closed_gc) < 3:
            if vassiliev:
                return 1, 0
            else:
                return 1, 1, 1

        alex_closed = alexander(closed_gc, -1)
        alex_closed = int(n.round(n.abs(alex_closed)))

        if vassiliev:
            v2_closed = vassiliev_degree_2(closed_gc)
            v2_closed = int(n.round(n.abs(v2_closed)))
            return alex_closed, v2_closed
        else:
            root_2 = n.exp(2 * n.pi * 1.j / 3.)
            root_2 = root_2.real.round(3) + root_2.imag.round(3)*1.j
            root_3 = n.exp(2 * n.pi * 1.j / 4.)
            root_3 = root_3.real.round(3) + root_3.imag.round(3)*1.j

            alex_closed_2 = alexander(closed_gc, root_2)
            alex_closed_2 = int(n.round(n.abs(alex_closed_2)))
            alex_closed_3 = alexander(closed_gc, root_3)
            alex_closed_3 = int(n.round(n.abs(alex_closed_3)))
            return alex_closed, alex_closed_2, alex_closed_3

    def closure_invariants(self, number_of_samples=10,
                              zero_centroid=False, vassiliev=False, **kwargs):
        '''
        '''
        if zero_centroid:
            self.zero_centroid()

        angles = get_rotation_angles(number_of_samples)

        full_invariants = []

        print_dist = int(max(1, 3000. / len(self.points)))
        for i, angs in enumerate(angles):
            if i % print_dist == 0:
                self._vprint('\ri = {} / {}'.format(i, len(angles)), False)
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            invariant = k.closure_invariant(vassiliev=vassiliev, **kwargs)
            if vassiliev:
                full_invariants.append([angs[0], angs[1], invariant[0], invariant[1]])
            else:
                full_invariants.append([angs[0], angs[1], invariant[0], invariant[1], invariant[2]])

        return n.array(full_invariants)

    def closure_invariants_fractions(self, number_of_samples=10, summary=False, vassiliev=False, **kwargs):
        '''
        '''
        polys = self.closure_invariants(number_of_samples, vassiliev=vassiliev, **kwargs)
        full_invariants = polys[:, 2:].astype(n.int)

        unique_invariants = [full_invariants[0]]
        for rows in full_invariants:
            unique_check = True
            for entries in unique_invariants:
                if (entries == rows).all():
                    unique_check = False
            if unique_check:
                unique_invariants.append(rows)

        fracs = []
        length = float(len(full_invariants))

        for uniques in unique_invariants:
            occurences = 0
            for alex in full_invariants:
                if (alex == uniques).all():
                    occurences += 1
            fracs.append((list(uniques), occurences / length))

        full_invariant_fracs = sorted(fracs, key=lambda j: j[1])

        if summary:
            knot_percent = 0.0

            for closure_info in full_invariant_fracs:
                if ((vassiliev and closure_info[0][0] != 1 and closure_info[0][1] != 0) or
                   closure_info[0] != [1, 1, 1]):

                    knot_percent += closure_info[1]

            return full_invariant_fracs, knot_percent
        else:
            return full_invariant_fracs

    def invariants_to_colours(self, theta=0, phi=0, virtual=True, jones=True):
        '''
        '''
        colours = {'2': (11/255., 247/255., 1., 1.0),  # cyan
                   '4': (54/255., 57/255., 1., 1.0),  # blue
                   '5': (2/255., 232/255., 102/255., 1.0),  # soft green
                   '[3, 4, 5]': (11/255., 247/255., 1., 1.0),  # cyan
                   '[3, 8, 5]': (140/255., 183/255., 245/255., 1.0),  # soft slate blue
                   '[7, 8, 9]': (61/255., 232/255., 11/255., 1.0),  # lime
                   '[3, 8, 9]': (157/255., 1., 213/255., 1.0),  # washed out cyan
                   '[7, 8, 19]': (146/255., 232/255., 143/255., 1.0),  # washed out green
                   '[1, 0]': (0.5, 0.5, 0.5, 1.0),  # grey
                   '[3, 1]': (1., 61/255., 0, 1.0),  # orange
                   '[5, 1]': (1., 0., 47/255., 1.0),  # fairly pink
                   '[5, 3]': (1., 209/255., 0., 1.0),  # yellow
                   '[7, 2]': (1., 102/255., 97/255., 1.0),  # washed out red
                   '[9, 2]': (232/255., 136/255., 89/255., 1.0),  # washed out orange
                   '[11, 1]': (236/255., 1., 97/255., 1.0),  # washed out sand
                   '[13, 1]': (232/255., 207/255., 89/255., 1.0)}  # washed out yellow/brown
        other_virtual = (0, 0, 0, 1)  # black
        other_classical = (1, 1, 1, 1)  # white
        if virtual:
            if jones:
                invariant = self.jones_for_distinguishing(theta=theta, phi=phi, vassiliev=True)  # vassiliev=True so that code does not need to change. good enough for all the knots we are distinguishing here anytway
            else:
                invariant = self.projection_invariant(theta=theta, phi=phi, vassiliev=True)
        else:
            invariant = self.closure_invariant(theta=theta, phi=phi, vassiliev=True)
            invariant = [int(invariant[0]), int(invariant[1])]

        if len(invariant) != 2:
            if invariant[1] == '-':
                invariant = int(float(invariant[0]))
            #elif int(float(invariant[2])) != 0:

            # elif invariant[0:3] != [0, 0, 0]:
            #     invariant = invariant[0:3]
            #     invariant = [int(x) for x in invariant]
            # else:
            #     invariant = invariant[3:5]
            #     invariant = [int(x) for x in invariant]

            elif invariant[3:] == [0, 0, 0, 0] or (invariant[3:5] != invariant[5:7]):
                invariant = invariant[0:3]
                invariant = [int(x) for x in invariant]
            else:
                invariant = invariant[3:5]
                invariant = [int(x) for x in invariant]

        if str(invariant) in colours:
            return colours[str(invariant)]
        elif type(invariant) == int:
            return other_virtual
        elif len(invariant) == 2:
            return other_classical
        else:
            return other_virtual

    def vassiliev_degree_2_average(self, samples=10, recalculate=False,
                                   **kwargs):
        '''Returns the average Vassliev degree 2 invariant calculated by
        averaging its combinatorial value over many different
        projection directions.

        Parameters
        ----------
        samples : int
            The number of directions to average over. Defaults to 10.
        recalculate : bool
            Whether to recalculate the writhe even if a cached result
            is available. Defaults to False.
        **kwargs :
            These are passed directly to :meth:`raw_crossings`.
        '''
        if (self._cached_v2 and
            samples in self._cached_v2 and
            not recalculate):
            return self._cached_v2[samples]

        from pyknot2.spacecurves.rotation import get_rotation_angles, rotate_to_top
        angles = get_rotation_angles(samples)

        v2s = []
        for theta, phi in angles:
            k = OpenKnot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(theta, phi))
            v2 = k._vassiliev_degree_2_projection()
            v2s.append(v2)

        result = n.average(v2s)
        self._cached_v2[samples] = result
        return result

    def _vassiliev_degree_2_projection(self):
        from pyknot2.spacecurves import Knot
        k = Knot(self.points, verbose=False)
        return k.vassiliev_degree_2(simplify=False,
                                    include_closure=False)


def gall_peters(theta, phi):
    '''
    Converts spherical coordinates to the Gall-Peters
    projection of the sphere, an area-preserving projection in
    the shape of a Rectangle.

    Parameters
    ----------
    theta : float
        The latitude, in radians.
    phi : float
        The longitude, in radians.
    '''
    theta -= n.pi / 2
    return (phi, 2 * n.sin(theta))

def knot_db_to_string(database_object):
    '''
    Takes output from from_invariants() and returns knot type as decimal.
    For example: <Knot 3_1> becomes 3.1 and <Knot K13n1496> becomes 13.1496
    '''
    from pyknot2.catalogue.identify import from_invariants
    db_strings = []
    for entries in database_object:
        db_string = str(entries)
        if db_string[6] == 'K':
            db_string = db_string[7:-1]
        else:
            db_string = db_string[6:-1]
            db_string = db_string.replace('_', '.')
        db_strings.append(db_string)
    return db_strings

