'''
OpenKnot
========

A class for working with the topology of open curves.
'''

from __future__ import print_function
import numpy as n
import sympy as sym

from pyknot2.spacecurves.spacecurve import SpaceCurve
from pyknot2.spacecurves.knot import Knot
from pyknot2.spacecurves.rotation import get_rotation_angles, rotate_to_top
from pyknot2.catalogue.identify import from_invariants
from pyknot2.catalogue.database import Knot as dbknot
from collections import Counter


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

    @property
    def points(self):
        return super(OpenKnot, self).points

    @points.setter
    def points(self, points):
        super(OpenKnot, self.__class__).points.fset(self, points)
        self._cached_alexanders = None

    def raw_crossings(self, mode='use_max_jump',
                      recalculate=False, use_python=False):
        '''
        Calls :meth:`pyknot2.spacecurves.spacecurve.SpaceCurve.raw_crossings`,
        but without including the closing line between the last
        and first points (i.e. setting include_closure=False).
        '''
        return super(OpenKnot, self).raw_crossings(mode=mode,
                                                   include_closure=False,
                                                   recalculate=recalculate,
                                                   use_python=use_python)

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

    def _plot_projections(self, number_of_samples):
        '''
        Plots the projection of the knot at each of the given
        number of samples, approximately evenly distributed on
        the sphere.

        This function is really just for testing.

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
                              zero_centroid=False):
        '''
        Returns a list of Alexander polynomials for the knot, closing
        on a sphere of the given radius, with the given number of sample
        points approximately evenly distributed on the sphere.

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
        if zero_centroid:
            self.zero_centroid()

        if self._cached_alexanders is not None:
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
            k = Knot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            points = k.points
            closure_point = points[-1] + points[0] / 2.
            closure_point[2] = radius
            k.points = n.vstack([points, closure_point])
            polys.append([angs[0], angs[1], k.alexander_polynomial()])

        self._cached_alexanders[
            (number_of_samples, cache_radius)] = n.array(polys)

        return n.array(polys)


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
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def _alexander_map_values(self, number_of_samples=10, **kwargs):
        polys = self.alexander_polynomials(
            number_of_samples=number_of_samples, **kwargs)

        from scipy.interpolate import griddata

        positions = []
        for i, row in enumerate(polys):
            positions.append(gall_peters(row[0], row[1]))
            '''positions.append((row[0]-n.pi/2, row[1]-n.pi))'''
        positions = n.array(positions)

        interpolation_points = n.mgrid[0:2 * n.pi:157j,
                               -2.:2.:100j]
        '''interpolation_points = n.mgrid[-n.pi:n.pi:157j,
                               -n.pi/2:n.pi/2:100j]'''

        values = griddata(positions, polys[:, 2],
                          tuple(interpolation_points),
                          method='nearest')

        return positions, values


    def plot_alexander_map(self, number_of_samples=10,
                           scatter_points=False,
                           mode='imshow', **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different Alexander polynomial in a different colour according
        to a closure on a far away point in this direction.
        '''

        positions, values = self._alexander_map_values(number_of_samples)

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

        self.plot()

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

    def virtual_check(self):
        '''
        Takes an open curve and checks (for the default projection) if its
        Gauss code corresponds to a virtual knot or not. Returns a Boolean of
        this information.

        Returns
        -------
        : virtual : bool
            True if the Gauss code corresponds to a virtual knot. False
            otherwise
        '''
        gauss_code = self.gauss_code()._gauss_code[0][:, 0]
        l = len(gauss_code)
        total_crossings = l / 2
        crossing_counter = 1
        virtual = False

        for crossing_number in self.gauss_code().crossing_numbers:
            occurences = n.where(gauss_code == crossing_number)[0]
            first_occurence = occurences[0]
            second_occurence = occurences[1]
            crossing_difference = second_occurence - first_occurence

            if crossing_difference % 2 == 0:
                return True
        return False


    def virtual_check_projections(self, number_of_samples=10,
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
            isvirtual = k.virtual_check()
            polys.append([angs[0], angs[1], isvirtual])

        return n.array(polys)

    def virtual_fractions(self, number_of_samples=10, **kwargs):
        '''Returns each of the virtual Booleans from
        self.virtual.check.projections, with the fraction of each type.
        '''
        polys = self.virtual_check_projections(
            number_of_samples=number_of_samples, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def _virtual_map_values(self, number_of_samples=10, **kwargs):
        polys = self.virtual_check_projections(
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

    def plot_virtual_map(self, number_of_samples=10,
                         scatter_points=False,
                         mode='imshow', **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different virtual Boolean in a different colour according
        to a projection in this direction.
        '''

        positions, values = self._virtual_map_values(number_of_samples)

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        if mode == 'imshow':
            cax = ax.imshow(values.T, cmap='jet', interpolation='none')
            fig.colorbar(cax)
        else:
            ax.contourf(values.T, cmap='jet',
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

    def plot_virtual_shell(self, number_of_samples=10,
                           zero_centroid=False,
                           sphere_radius_factor=2.,
                           opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured according to whether or not the projection from this
        point corresponds to a virtual knot or not.

        Parameters are all passed to
        :meth:`OpenKnot.virtual_check_projections`, except opacity and kwargs
        which are given to mayavi.mesh, and sphere_radius_factor which gives
        the radius of the enclosing sphere in terms of the maximum Cartesian
        distance of any point in the line from the origin.
        '''

        self.plot()

        positions, values = self._virtual_map_values(
            number_of_samples, zero_centroid=False)

        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
        phis = n.linspace(0, 2 * n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor * n.max(self.points)
        zs = r * n.cos(thetas)
        xs = r * n.sin(thetas) * n.cos(phis)
        ys = r * n.sin(thetas) * n.sin(phis)

        import mayavi.mlab as may

        may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)


    def self_link(self):
        '''
        Takes an open curve, finds its Gauss code (for the default projection)
        and calculates its self linking number, J(K). See Kauffman 2004 for
        more information.

        Returns
        -------
        : self_link_counter : int
            The self linking number of the open curve
        '''

        gausscode = self.gauss_code()._gauss_code
        l = len(gausscode[0][:, 0])
        total_crossings = l / 2
        crossing_number = 1
        self_link_counter = 0

        for crossing_number in self.gauss_code().crossing_numbers:
            occurrences = n.where(gausscode[0][:, 0] == crossing_number)[0]
            first_occurrence = occurrences[0]
            second_occurrence = occurrences[1]
            crossing_difference = second_occurrence - first_occurrence

            if (crossing_difference % 2 == 0):
                self_link_counter += gausscode[0][occurrences[0], 2]

        return self_link_counter/2


    def self_link_projections(self, number_of_samples=10,
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
            self_link = k.self_link()
            polys.append([angs[0], angs[1], self_link])

        return n.array(polys)


    def self_link_fractions(self, number_of_samples=10, **kwargs):
        '''Returns each of the self linking numbers from
        self.virtual.self_link.projections, with the fraction of each type.
        '''
        polys = self.self_link_projections(
            number_of_samples=number_of_samples, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]


    def _self_link_map_values(self, number_of_samples=10, **kwargs):
        polys = self.self_link_projections(
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


    def plot_self_link_map(self, number_of_samples=10,
                           scatter_points=False,
                           mode='imshow', **kwargs):
        '''
        Creates (and returns) a projective diagram showing each
        different self linking number in a different colour according
        to a projection in this direction.
        '''

        positions, values = self._self_link_map_values(number_of_samples)

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        if mode == 'imshow':
            cax = ax.imshow(values.T, cmap='jet', interpolation='none')
            fig.colorbar(cax)
        else:
            ax.contourf(values.T, cmap='jet',
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


    def plot_self_link_shell(self, number_of_samples=10,
                             zero_centroid=False,
                             sphere_radius_factor=2.,
                             opacity=0.3, **kwargs):
        '''
        Plots the curve in 3d via self.plot(), along with a translucent
        sphere coloured by the self linking number obtained by projecting from
        this point.

        Parameters are all passed to
        :meth:`OpenKnot.virtual_check_projections`, except opacity and kwargs
        which are given to mayavi.mesh, and sphere_radius_factor which gives
        the radius of the enclosing sphere in terms of the maximum Cartesian
        distance of any point in the line from the origin.
        '''

        self.plot()

        positions, values = self._self_link_map_values(
            number_of_samples, zero_centroid=False)

        thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
        phis = n.linspace(0, 2 * n.pi, 157)

        thetas, phis = n.meshgrid(thetas, phis)

        r = sphere_radius_factor * n.max(self.points)
        zs = r * n.cos(thetas)
        xs = r * n.sin(thetas) * n.cos(phis)
        ys = r * n.sin(thetas) * n.sin(phis)

        import mayavi.mlab as may

        may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)


    def alexander_polynomials_multiroots(self, number_of_samples=10, radius=None,
                                         zero_centroid=False):
        '''
        Returns a list of Alexander polynomials for the knot, closing
        on a sphere of the given radius, with the given number of sample
        points approximately evenly distributed on the sphere. The
        Alexander polynomials are found at three different roots (2, 3
        and 4) and a the knot types corresponding to these roots are
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
        if zero_centroid:
            self.zero_centroid()

        if self._cached_alexanders is not None:
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
            k_gauss_code.simplify(verbose = False)
            max_crossings = len(k_gauss_code)
            if max_crossings > 17:
                max_crossings = 18
            knot_type = from_invariants(determinant=root_at_two, alex_imag_3=root_at_three,
                                        alex_imag_4=root_at_four,
                                        other=[dbknot.min_crossings <= max_crossings])
            knot_type = knot_db_to_string(knot_type)
            polys.append([angs[0], angs[1], root_at_two, root_at_three, root_at_four, max_crossings,
                          knot_type])

            # self._cached_alexanders[
            #      (number_of_samples, cache_radius)] = n.array(polys)

        return polys

    def multiroots_fractions(self, number_of_samples=10, **kwargs):

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

    def generalised_alexander(self):
        '''
        Returns the generalised Alexander polynomial for the default projection
        of the open knot
        '''

        gauss_code_crossings = self.gauss_code()._gauss_code[0][:, 0]
        gauss_code_over_under = self.gauss_code()._gauss_code[0][:,1]
        gauss_code_orientations = self.gauss_code()._gauss_code[0][:,2]
        x = sym.var('x')
        y = sym.var('y')
        m_plus = sym.Matrix([[1 - x, -y], [-x * y**-1, 0]])
        m_minus = sym.Matrix([[0, -x**-1 * y], [-y**-1, 1 - x**-1]])
        num_crossings = len(self.gauss_code())
        matrix = sym.zeros(2 * num_crossings, 2 *num_crossings)
        permutation_matrix = sym.zeros(2 * num_crossings, 2 * num_crossings)

        arc_labels = [0]*len(gauss_code_crossings)
        for i in range(len(gauss_code_crossings)):
            arc_labels[i] = [0] * 4

        counter = 0
        for crossing_number in self.gauss_code().crossing_numbers:
            occurrences = n.where(gauss_code_crossings == crossing_number)[0]
            if gauss_code_orientations[occurrences[0]] == 1:
                m = m_plus
            else:
                m = m_minus
            for i in [0,1]:
                for j in [0,1]:
                    matrix[counter*2 + i,counter*2 + j] = m[i,j]
            counter += 1

        for i in range(len(gauss_code_crossings)):
            arc_labels[i][0] = gauss_code_crossings[i]
            arc_labels[i][1] = (gauss_code_orientations[i] *
                                gauss_code_over_under[i]) #-1 = r, +1 = l
            arc_labels[i-1][2] = gauss_code_crossings[i]
            arc_labels[i-1][3] = (-1 * gauss_code_orientations[i] *
                                  gauss_code_over_under[i]) #-1 = r, +1 = l

        counter = 1
        for crossing_number in self.gauss_code().crossing_numbers:
            for i in range(len(gauss_code_crossings)):
                if arc_labels[i][0] == crossing_number:
                    arc_labels[i][0] = counter
                if arc_labels[i][2] == crossing_number:
                    arc_labels[i][2] = counter
            counter += 1

        for i in range(len(arc_labels)):
            if arc_labels[i][1] < 0:
                if arc_labels[i][3] < 0 :
                    #bottom right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-1] = 1
                else:
                    #bottom left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-1] = 1
            else:
                if arc_labels[i][3] < 0:
                    #upper right
                    permutation_matrix[arc_labels[i][2]*2-1, arc_labels[i][0]*2-2] = 1
                else:
                    #upper left
                    permutation_matrix[arc_labels[i][2]*2-2, arc_labels[i][0]*2-2] = 1

        writhe = sum(gauss_code_orientations)/2

        return (-1)**writhe * ((matrix - permutation_matrix).det())

    def projection_invariant(self, **kwargs):
        '''
        First checks if the projection of an open curve is virtual or classical. If virtual,
        a virtual knot invariant is calculated. Otherwise a classical invariant is calculated.
        '''
        is_virtual = self.virtual_check()
        if is_virtual == True:
            return(['v', self.self_link()])
        else:
            from pyknot2.invariants import alexander
            return(['c'],alexander(self.gauss_code(), variable=-1, quadrant='lr',
                                   mode='python', simplify=False))

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

def mollweide(phi, lambda_):
    '''

    Converts spherical coordinates to the Mollweide
    projection of the sphere, an area-preserving projection in
    the shape of an ellipse.

    Parameters
    ----------
    phi : float
        The latitude, in radians.
    lambda_ : float
        The longitude, in radians.
    '''

    if phi == n.pi/2 or phi == -n.pi/2:
        theta_m = phi
    else:
        theta_n = phi
        theta_m = 0
        while abs(theta_m - theta_n) > 0.000001:
            theta_n = theta_m
            theta_m = (theta_n - (2*theta_n + n.sin(2*theta_n) - n.pi * n.sin(phi)) /
                       (2 + 2*n.cos(2*theta_n)))
    x = ((2 * n.sqrt(2)) / (n.pi)) * lambda_ * n.cos(theta_m)
    y = n.sqrt(2) * n.sin(theta_m)
    return(x,y)


def knot_db_to_string(database_object):
    '''
    Takes output from from_invariants() and returns knot type as decimal.
    For example: <Knot 3_1> becomes 3.1 and <Knot K13n1496> becomes 13.1496
    '''
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


