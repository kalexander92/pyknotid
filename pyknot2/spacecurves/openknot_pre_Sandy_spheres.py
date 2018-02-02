'''
OpenKnot
========

A class for working with the topology of open curves.
'''

from __future__ import print_function
import numpy as n

from pyknot2.spacecurves.spacecurve import SpaceCurve
from pyknot2.spacecurves.knot import Knot
from pyknot2.spacecurves.rotation import get_rotation_angles, rotate_to_top
from pyknot2.catalogue.identify import from_invariants
from pyknot2.catalogue.database import Knot as dbknot
from collections import Counter

from pyknot2.invariants import alexander, generalised_alexander, vassiliev_degree_2
from pyknot2.representations.gausscode import GaussCode


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
                      recalculate=False, try_cython=False):
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
            closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-1)) |
                                    ((cs[:, 1] > len(self.points)-1)))
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
                              optimise_closure=True):
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
            k = Knot(self.points, verbose=False)
            k._apply_matrix(rotate_to_top(*angs))
            if zero_centroid:
                k.zero_centroid()
            if optimise_closure:
                cs = k.raw_crossings()
                if len(cs) > 0:
                    closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-1) & (cs[:, 2] < 0.)) |
                                            ((cs[:, 1] > len(self.points)-1) & (cs[:, 2] > 0.)))
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

                polys.append([angs[0], angs[1], k.alexander_polynomial(**kwargs)])

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
            number_of_samples, radius=None, zero_centroid=False)

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
        '''Returns each of the virtual booleans from
        self.virtual.check.projections, with the fraction of each type.
        '''
        polys = self.virtual_checks(
            number_of_samples=number_of_samples, **kwargs)
        alexs = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(alexs))
        for alex in n.unique(alexs):
            fracs.append((alex, n.sum(alexs == alex) / length))

        return sorted(fracs, key=lambda j: j[1])
        
    def _virtual_map_values(self, number_of_samples=10, **kwargs):
        polys = self.virtual_checks(
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
        :meth:`OpenKnot.virtual_checks`, except opacity and kwargs 
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

    def self_linking(self):
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
        crossing_number = 1
        self_link_counter = 0

        for crossing_number in self.gauss_code().crossing_numbers:
            occurrences = n.where(gausscode[0][:, 0] == crossing_number)[0]
            first_occurrence = occurrences[0]
            second_occurrence = occurrences[1]
            crossing_difference = second_occurrence - first_occurrence

            if (crossing_difference % 2 == 0):
                self_link_counter += gausscode[0][occurrences[0], 2]

        return self_link_counter

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
            for i in [0, 1]:
                for j in [0, 1]:
                    matrix[counter*2 + i,counter*2 + j] = m[i,j]
            counter += 1

        for i in range(len(gauss_code_crossings)):
            arc_labels[i][0] = gauss_code_crossings[i]
            arc_labels[i][1] = (gauss_code_orientations[i] *
                                gauss_code_over_under[i])  # -1 = r, +1 = l
            arc_labels[i-1][2] = gauss_code_crossings[i]
            arc_labels[i-1][3] = (-1 * gauss_code_orientations[i] *
                                  gauss_code_over_under[i])  # -1 = r, +1 = l

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
        for crossing_number in self.gauss_code().crossing_numbers:
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
        for crossing_number in self.gauss_code().crossing_numbers:
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

    def jones_polynomial(self, root=-1, simplify=True):
        '''
        Calculates the jones polynomial for a single projection of the open curve
        using the given root
        '''
        import pyknot2.representations.planardiagram as pdiag
        if simplify:
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
            jones = k.jones_polynomial(root, simplify)
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

    def projection_invariant(self):
        '''
        First calculates the generalised Alexander polynomial of the curve. If 0,
        then calculates the classical Alexander determinant and vassiliev invariant
        of degree 2.

        Returned is a list of invariants in the following order:
        [gen_alex(2,3), gen_alex(2,4), gen_alex(3,4), alex(-1), vassiliev_2]
        Returns 0 for classical invariants if projection is virtual
        '''
        k = Knot(self.points, verbose=False)
        cs = k.raw_crossings()
        if len(cs) > 0:
            closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-1)
                                    & (cs[:, 2] < 0.)) |
                                    ((cs[:, 1] > len(self.points)-1)
                                    & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            for index in indices:
                cs[index, 2:] *= -1
        closed_gc = GaussCode(cs, verbose=self.verbose)

        # Remove closing crossings to calculate self linking
        xs, ys = n.where(cs[:, :2] > len(k.points) - 1)
        keeps = n.ones(len(cs), dtype=n.bool)
        for x in xs:
            keeps[x] = False
        cs = cs[keeps]
        open_gc = GaussCode(cs)

        gen_alex_1 = generalised_alexander(open_gc, root_s=2, root_t=3)
        if gen_alex_1 != 0:
            gen_alex_2 = generalised_alexander(open_gc, root_s=2, root_t=4)
            gen_alex_3 = generalised_alexander(open_gc, root_s=3, root_t=4)
            return [gen_alex_1, gen_alex_2, gen_alex_3, 0, 0, 0, 0]
        else:
            #closed_gc.simplify()
            alex_closed = alexander(closed_gc, -1)
            alex_open = alexander(open_gc, -1)
            alex_closed = int(n.round(n.abs(alex_closed)))
            alex_open = int(n.round(n.abs(alex_open)))

            v2_closed = vassiliev_degree_2(closed_gc)
            v2_open = vassiliev_degree_2(open_gc)
            v2_closed = int(n.round(n.abs(v2_closed)))
            v2_open = int(n.round(n.abs(v2_open)))
            return [0, 0, 0, alex_closed, v2_closed, alex_open, v2_open]

    def projection_invariants(self, number_of_samples=10,
                              zero_centroid=False):
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
            invariants = k.projection_invariant()
            full_invariants.append([angs[0], angs[1], invariants[0], invariants[1],
                                    invariants[2], invariants[3], invariants[4], invariants[5],
                                    invariants[6]])

        return n.array(full_invariants)

    def projection_invariants_fractions(self, number_of_samples=10, summary=False):
        '''
        Returns each of the projection invariants from
        self.virtual.projections_invariants, with the fraction of each type. If
        summary=True also returns [unknot %, classical %, virtual %],
        [most common classical invariant, %] and
        [most common virtual invariant, %]
        '''
        polys = self.projection_invariants(number_of_samples)
        full_invariants = polys[:, 2:].astype(n.int)
        closure_invariants = polys[:, 2:7].astype(n.int)

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
            unique_invariants = [closure_invariants[0]]
            for rows in closure_invariants:
                unique_check = True
                for entries in unique_invariants:
                    if (entries == rows).all():
                        unique_check = False
                if unique_check:
                    unique_invariants.append(rows)

            fracs = []
            length = float(len(closure_invariants))

            for uniques in unique_invariants:
                occurences = 0
                for alex in closure_invariants:
                    if (alex == uniques).all():
                        occurences += 1
                fracs.append((list(uniques), occurences / length))

            closure_invariant_fracs = sorted(fracs, key=lambda j: j[1])

            classical_number = 0.0
            virtual_number = 0.0
            unknot_number = 0.0
            highest_classical = 0
            classical_mode = (None, None)
            highest_virtual = 0
            virtual_mode = (None, None, None)

            for projection_info in closure_invariant_fracs:
                if projection_info[0][3] == 1 and projection_info[0][4] == 0:
                    unknot_number += projection_info[1]
                elif projection_info[0][0] == 0:
                    classical_number += projection_info[1]
                    if projection_info[1] >= highest_classical:
                        classical_mode = projection_info[0]
                        highest_classical = projection_info[1]
                else:
                    virtual_number += projection_info[1]
                    if projection_info[1] >= highest_virtual:
                        virtual_mode = projection_info[0]
                        highest_virtual = projection_info[1]
            return (full_invariant_fracs, [unknot_number, classical_number, virtual_number],
                    [classical_mode, highest_classical], [virtual_mode, highest_virtual])
        else:
            return full_invariant_fracs

    def jones_for_distinguishing(self, root=-1, simplify=True):
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

        alexander_check = self.projection_invariant()
        if alexander_check == [3, 4, 5, 0, 0, 0, 0]:
            return [self.jones_polynomial(root, simplify), '-', '-', '-', '-', '-', '-']
        else:
            return alexander_check

    def jones_for_distinguishings(self, number_of_samples=10, root=-1, simplify=True, **kwargs):
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
            jones = k.jones_for_distinguishing(root, simplify)
            polys.append([angs[0], angs[1], jones[0], jones[1], jones[2], jones[3], jones[4], jones[5], jones[6]])

        return n.array(polys)

    def jones_for_distinguishings_fractions(self, number_of_samples=10, root=-1, simplify=True, **kwargs):
        '''
        Returns each of the jones_for_distinguishing with the fraction of each.
        '''
        '''
        polys = self.jones_for_distinguishings(
            number_of_samples=number_of_samples, root=root, simplify=True, **kwargs)
        joneses = n.round(polys[:, 2]).astype(n.int)

        fracs = []
        length = float(len(joneses))
        for jones in n.unique(joneses):
            fracs.append((jones, n.sum(joneses == jones) / length))
        # fracs = n.array(fracs)

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]
        '''

        polys = self.jones_for_distinguishings(
            number_of_samples=number_of_samples, root=root, simplify=simplify, **kwargs)
        joneses = polys[:, 2:9]

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

        return sorted(fracs, key=lambda j: j[1])
        #return fracs[n.argsort(fracs[:, 1])]

    def closure_invariant(self):
        '''
        '''
        k = Knot(self.points, verbose=False)
        cs = k.raw_crossings()
        if len(cs) > 0:
            closure_cs = n.argwhere(((cs[:, 0] > len(self.points)-1)
                                    & (cs[:, 2] < 0.)) |
                                    ((cs[:, 1] > len(self.points)-1)
                                    & (cs[:, 2] > 0.)))
            indices = closure_cs.flatten()
            for index in indices:
                cs[index, 2:] *= -1
        closed_gc = GaussCode(cs, verbose=self.verbose)

        closed_gc.simplify()
        alex_closed = alexander(closed_gc, -1)
        alex_closed = int(n.round(n.abs(alex_closed)))

        v2_closed = vassiliev_degree_2(closed_gc)
        v2_closed = int(n.round(n.abs(v2_closed)))
        return alex_closed, v2_closed

    def closure_invariants(self, number_of_samples=10,
                              zero_centroid=False):
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
            alex, v2 = k.closure_invariant()
            full_invariants.append([angs[0], angs[1], alex, v2])

        return n.array(full_invariants)

    def closure_invariants_fractions(self, number_of_samples=10, summary=False):
        '''
        '''
        polys = self.closure_invariants(number_of_samples)
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
                if closure_info[0][0] != 1 and closure_info[0][1] != 0:

                    knot_percent += closure_info[1]

            return full_invariant_fracs, knot_percent
        else:
            return full_invariant_fracs

    def plot_geb_figure(self, spacings=(1, 1, 1), panel=False):
        '''
        '''

        self.zero_centroid()

        k_min = n.amin(self.points, axis=0)
        k_max = n.amax(self.points, axis=0)

        ranges = [abs(k_max[x] - k_min[x]) for x in range(3)]

        translation_vector = [abs(k_min[x]) + n.min(ranges)*spacings[x] for x in range(3)]

        self.translate(translation_vector)

        k_x = OpenKnot([[0, x[1], x[2]] for x in self.points])
        k_y = OpenKnot([[x[0], 0, x[2]] for x in self.points])
        k_z = OpenKnot([[x[0], x[1], 0] for x in self.points])

        if not panel:
            self.plot(mode='vispy', zero_centroid=False, colour=(0/255.0, 89/255.0, 255/255.0))

            import pyknot2.visualise as vis
            vis.plot_vispy_cube(dimensions=n.average(k_x.points, axis=0), clf=False, colour=(1, 0.9, 0.8, 1),
                                edge_colour='black', offset=n.average(k_x.points, axis=0))
            vis.plot_vispy_cube(dimensions=n.average(k_y.points, axis=0), clf=False, colour=(0.8, 1, 0.8, 1),
                                edge_colour='black', offset=n.average(k_y.points, axis=0))
            vis.plot_vispy_cube(dimensions=n.average(k_z.points, axis=0), clf=False, colour=(0.8, 0.9, 1, 1),
                                edge_colour='black', offset=n.average(k_z.points, axis=0))

            k_x.plot(clf=False, zero_centroid=False, colour='black')
            k_y.plot(clf=False, zero_centroid=False, colour='black')
            k_z.plot(clf=False, zero_centroid=False, colour='black')

            vis.camera_position(fov=10, distance=20*n.max(ranges), elevation=30,
                                azimuth=145, centre=n.average(self.points, axis=0))
        if panel == 'x':
            import pyknot2.visualise as vis
            vis.plot_vispy_cube(dimensions=n.average(k_x.points, axis=0), colour=(1, 0.9, 0.8, 1),
                                edge_colour='black', offset=n.average(k_x.points, axis=0))
            k_x.plot(clf=False, zero_centroid=False, colour='black')
            vis.camera_position(fov=10, distance=20*n.max(ranges), elevation=0,
                                azimuth=90, centre=n.average(self.points, axis=0))

        if panel == 'y':
            import pyknot2.visualise as vis
            vis.plot_vispy_cube(dimensions=n.average(k_y.points, axis=0), colour=(0.8, 1, 0.8, 1),
                                edge_colour='black', offset=n.average(k_y.points, axis=0))
            k_y.plot(clf=False, zero_centroid=False, colour='black')
            vis.camera_position(fov=10, distance=20*n.max(ranges), elevation=0,
                                azimuth=180, centre=n.average(self.points, axis=0))

        if panel == 'z':
            import pyknot2.visualise as vis
            vis.plot_vispy_cube(dimensions=n.average(k_z.points, axis=0), colour=(0.8, 0.9, 1, 1),
                                edge_colour='black', offset=n.average(k_z.points, axis=0))
            k_z.plot(clf=False, zero_centroid=False, colour='black')
            vis.camera_position(fov=10, distance=20*n.max(ranges), elevation=90,
                                azimuth=90, centre=n.average(self.points, axis=0))

    def plot_geb_figure_panels(self, spacings=(1, 1, 1)):
        '''
        '''
        import pyknot2.visualise as vis
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg

        self.plot_geb_figure(spacings=spacings)
        vis.vispy_save_png('geb_full.png')
        self.plot_geb_figure(spacings=spacings, panel='x')
        vis.vispy_save_png('geb_x.png')
        self.plot_geb_figure(spacings=spacings, panel='y')
        vis.vispy_save_png('geb_y.png')
        self.plot_geb_figure(spacings=spacings, panel='z')
        vis.vispy_save_png('geb_z.png')

        img_full = mpimg.imread('geb_full.png')
        x_full = mpimg.imread('geb_x.png')
        y_full = mpimg.imread('geb_y.png')
        z_full = mpimg.imread('geb_z.png')

        ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=3)
        ax1.imshow(_crop_whitespace(img_full))
        ax1.axes.get_xaxis().set_visible(False)
        ax1.axes.get_yaxis().set_visible(False)

        ax2 = plt.subplot2grid((2, 3), (1, 0))
        ax2.imshow(_crop_whitespace(x_full))
        ax2.axes.get_xaxis().set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)

        ax3 = plt.subplot2grid((2, 3), (1, 1))
        ax3.imshow(_crop_whitespace(y_full))
        ax3.axes.get_xaxis().set_visible(False)
        ax3.axes.get_yaxis().set_visible(False)

        ax4 = plt.subplot2grid((2, 3), (1, 2))
        ax4.imshow(_crop_whitespace(z_full))
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)

        plt.tight_layout()


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

def _crop_whitespace(img):
    '''
    For cropping white space around panel figures in plot_geb_figure_panels
    '''
    x_len = len(img[0])
    y_len = len(img)
    #starting at opposite ends and moving borders past each other
    x_min = x_len
    y_min = y_len
    x_max = 0
    y_max = 0
    white = n.array([1.0, 1.0, 1.0, 1.0])

    for y in range(y_len):
        for x in range(x_len):
            if (img[y][x] != white).any():
                if x < x_min:
                    x_min = x
                elif x > x_max:
                    x_max = x
                if y < y_min:
                    y_min = y
                elif y > y_max:
                    y_max = y
    img = [x[x_min:x_max] for x in img]
    return img[y_min:y_max]






