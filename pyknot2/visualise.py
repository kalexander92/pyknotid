'''
Visualise
=========

This module contains functions for plotting knots, supporting
different toolkits and types of plot.
'''

import numpy as n
from colorsys import hsv_to_rgb
from pyknot2.utils import ensure_shape_tuple, vprint
import random
import signal

vispy_canvas = None

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

def camera_position(fov=30, distance=20, 
                    elevation=30, azimuth=45,
                    centre=[0, 0, 0]):
    ensure_vispy_canvas()
    canvas = vispy_canvas
    from vispy import app, scene, color

    canvas.view.camera = scene.TurntableCamera(fov=fov, distance=distance,
                                               elevation=elevation, azimuth=azimuth,
                                               center=centre)
    canvas.show()

def plot_vispy_cube(dimensions=(5, 2, 10), clf=True, 
                    colour=(0, 1, 0, 1), edge_colour=(0, 0, 0, 1),
                    offset = (0, 0, 0)):
    ensure_vispy_canvas()
    if clf:
        clear_vispy_canvas()
    canvas = vispy_canvas

    from vispy import scene, color

    c = scene.visuals.Cube(dimensions,
                           color=colour,
                           edge_color=edge_colour)
    c.transform = scene.transforms.MatrixTransform()
    #c.transform = scene.transforms.AffineTransform()
    c.transform.translate(offset)

    mesh = c.mesh
    md = mesh._meshdata
    vertices = md.get_vertices()
    faces = md.get_faces()

    colours = [colour for x in range(len(vertices))]
    md.set_vertex_colors(colours)
    colours = [edge_colour for x in range(len(faces))]
    md.set_face_colors(colours)


    ensure_vispy_canvas()
    canvas.view.camera = scene.ArcballCamera(fov=30, distance=20)
    canvas.view.add(c)
    canvas.show()

def plot_vispy_sphere(radius=1., rows=100, cols=100,
                      opacity=0.3,
                      method='latitude',
                      colour=(0, 0.5, 0.5, 0.3),
                      edge_color=None,
                      **kwargs):
    '''
    '''

    from vispy.scene import Sphere, ArcballCamera

    colour = (colour[0], colour[1], colour[2], opacity)

    s = Sphere(rows=rows, cols=cols, method=method,
               edge_color=edge_color, color=colour,
               radius=radius)

    mesh = s.mesh
    md = mesh._meshdata
    vertices = md.get_vertices()

    colours = [colour for x in range(len(vertices))]

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    vispy_canvas.view.camera = ArcballCamera(fov=30)
    camera_position(distance=int(radius) * 5)
    vispy_canvas.view.add(s)
    vispy_canvas.show()

def plot_line(points, mode='auto', clf=True, **kwargs):
    '''
    Plots the given line, using the toolkit given by mode.

    kwargs are passed to the toolkit specific function, except for:

    Parameters
    ----------
    points : ndarray
        The nx3 array to plot.
    mode : str
        The toolkit to draw with. Defaults to 'auto', which will
        automatically pick the first available toolkit from
        ['mayavi', 'matplotlib', 'vispy'], or raise an exception
        if none can be imported.
    clf : bool
        Whether the existing figure should be cleared
        before drawing the new one.
    '''
    if mode == 'auto':
        try:
            import vispy
            mode = 'vispy'
        except ImportError:
            pass
    if mode == 'auto':
        try:
            import mayavi.mlab as may
            mode = 'mayavi'
        except ImportError:
            pass
    if mode == 'auto':
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            mode = 'matplotlib'
        except (ImportError, ValueError):
            pass
    if mode == 'auto':
        raise ImportError('Couldn\'t import any of mayavi, vispy, '
                          'or matplotlib\'s 3d axes.')
            
    if mode == 'mayavi':
        plot_line_mayavi(points, clf=clf, **kwargs)
    elif mode == 'vispy':
        plot_line_vispy(points, clf=clf, **kwargs)
    elif mode == 'matplotlib':
        plot_line_matplotlib(points, clf=clf, **kwargs)
    else:
        raise Exception('invalid toolkit/mode')

    
def plot_line_mayavi(points, clf=True, tube_radius=1., colormap='hsv',
                     closed=True,
                     zero_centroid=False,
                     mus=None,
                     **kwargs):
    import mayavi.mlab as may
    if clf:
        may.clf()
    if mus is None:
        mus = n.linspace(0, 1, len(points))
    may.plot3d(points[:, 0], points[:, 1], points[:, 2], mus,
               colormap=colormap, tube_radius=tube_radius, **kwargs)

def plot_line_matplotlib(points, **kwargs):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(points[:, 0], points[:, 1], points[:, 2])
    fig.show()

def ensure_vispy_canvas():
    global vispy_canvas
    if vispy_canvas is None:
        from vispy import app, scene
        canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
        canvas.unfreeze()
        canvas.view = canvas.central_widget.add_view()
        vispy_canvas = canvas
    # if not vispy_canvas.central_widget.children:
    #     vispy_canvas.view = vispy_canvas.central_widget.add_view()
        

def clear_vispy_canvas():
    global vispy_canvas
    if vispy_canvas is None:
        return
    vispy_canvas.unfreeze()
    vispy_canvas.central_widget.remove_widget(vispy_canvas.view)
    vispy_canvas.view = vispy_canvas.central_widget.add_view()

def vispy_rotate(elevation=0, max_angle=360):
    global vispy_canvas
    cam = vispy_canvas.view.camera
    from vispy.scene import TurntableCamera
    from time import sleep
    vispy_canvas.view.camera = TurntableCamera(
        fov=cam.fov, scale_factor=cam.scale_factor,
        azimuth=0, elevation=elevation)
    try:
        for i in range(max_angle):
            vispy_canvas.view.camera.azimuth = i
            vispy_canvas.update()
            vispy_canvas.events.draw()
            vispy_canvas.swap_buffers()
            sleep(1 / 60.)
    except KeyboardInterrupt:
        pass
    vispy_canvas.view.camera = cam


def plot_line_vispy(points, clf=True, tube_radius=1.,
                    colour=None, zero_centroid=True,
                    closed=False, mus=None,
                    tube_points=8, **kwargs):
    # Add an extra point to fix tube drawing bug
    last_tangent = points[-1] - points[-2]
    points = n.vstack([points, points[-1] + 0.0001 * last_tangent])

    ensure_vispy_canvas()
    if clf:
        clear_vispy_canvas()
    canvas = vispy_canvas
    from vispy import app, scene, color

    if colour is None:
        from colorsys import hsv_to_rgb
        if closed:
            colours = n.linspace(0, 1, len(points))
        else:
            colours = n.linspace(0, 0.667, len(points))
        colours = n.array([hsv_to_rgb(c, 1, 1) for c in colours])
    else:
        colours = color.ColorArray(colour)

    if mus is not None:
        colours = n.array([hsv_to_rgb(c, 1, 1) for c in mus])

    black_ends = kwargs.get('black_ends', False)
    if black_ends:
        colours[0] = [0, 0, 0]
        colours[-1] = [0, 0, 0]
        colours[-2] = [0, 0, 0]
        colours[-3] = [0, 0, 0]

    l = scene.visuals.Tube(points, color=colours,
                           shading='smooth',
                           radius=tube_radius,
                           tube_points=tube_points,
                           closed=closed)

    canvas.view.add(l)
    # canvas.view.camera = 'arcball'
    canvas.view.camera = scene.ArcballCamera(fov=30, distance=7.5*n.max(
        n.abs(points)))
    #canvas.view.camera = scene.TurntableCamera(fov=30)
    if zero_centroid:
        l.transform = scene.transforms.MatrixTransform()
        # l.transform = scene.transforms.AffineTransform()
        l.transform.translate(-1*n.average(points, axis=0))

    canvas.show()
    # import ipdb
    # ipdb.set_trace()
    return canvas

def plot_lines_vispy(lines, clf=True, tube_radius=1.,
                     colours=None, zero_centroid=True, tube_points=8,
                     closed=False,
                     **kwargs):
    ensure_vispy_canvas()
    if clf:
        clear_vispy_canvas()
    canvas = vispy_canvas
    from vispy import app, scene, color

    if not isinstance(tube_radius, list):
        tube_radius = [tube_radius for _ in range(len(lines))]

    if colours is None:
        colours = ['purple' for line in lines]

    tubes = []
    for colour, points, radius in zip(colours, lines, tube_radius):

        l = scene.visuals.Tube(points, color=colour,
                               shading='smooth',
                               radius=radius,
                               closed=closed,
                               tube_points=tube_points)
        tubes.append(l)

    from visualcollection import MeshCollection
    collection = MeshCollection(tubes)
    canvas.view.add(collection)
    canvas.view.camera = 'arcball'
    canvas.view.camera.fov = 30
    # canvas.view.camera = scene.TurntableCamera(
    #     fov=90, up='z', distance=1.2*n.max(n.max(
    #         points, axis=0)))

    if zero_centroid:
        l.transform = scene.transforms.MatrixTransform()
        # l.transform = scene.transforms.AffineTransform()
        l.transform.translate(-1*n.average(points, axis=0))

    canvas.show()
    return canvas

def plot_line_with_crossings(points, raw_crossings, crossing_width=2,
                             clf=True, tube_radius=1.,
                             colour=None, zero_centroid=True,
                             closed=False, mus=None,
                             tube_points=8, **kwargs):
    # Add an extra point to fix tube drawing bug
    last_tangent = points[-1] - points[-2]
    points = n.vstack([points, points[-1] + 0.0001 * last_tangent])

    ensure_vispy_canvas()
    if clf:
        clear_vispy_canvas()
    canvas = vispy_canvas
    from vispy import scene, color

    if colour is None:
        from colorsys import hsv_to_rgb
        colours = n.linspace(0, 1, len(points))
        colours = n.array([hsv_to_rgb(c, 1, 1) for c in colours])
    else:
        colours = color.ColorArray(colour)

    if mus is not None:
        colours = n.array([hsv_to_rgb(c, 1, 1) for c in mus])

    l = scene.visuals.Tube(points, color=colours,
                           shading='smooth',
                           radius=tube_radius,
                           tube_points=tube_points,
                           closed=closed)

    undercrossings = []
    for crossing in raw_crossings:
        if crossing[2] == -1:
            undercrossings.append(int(crossing[0]))

    trans_cs = l.mesh_data.get_vertex_colors()
    for crossing in undercrossings:
        for i in range(tube_points):
            for j in range(crossing_width):
                trans_cs[(crossing + 1 - j) * tube_points + i] = n.array([0., 0., 0., 0.])
                trans_cs[(crossing + j) * tube_points + i] = n.array([0., 0., 0., 0.])

    l.mesh_data.set_vertex_colors(trans_cs)

    canvas.view.add(l)
    canvas.view.camera = scene.ArcballCamera(fov=30, distance=7.5*n.max(
        n.abs(points)))
    if zero_centroid:
        l.transform = scene.transforms.MatrixTransform()
        l.transform.translate(-1*n.average(points, axis=0))

    canvas.show()
    return canvas


def plot_projection(points, crossings=None, mark_start=False,
                    fig_ax=None, show=True):
    '''
    Plot the 2d projection of the given points, with optional
    markers for where the crossings are.

    Parameters
    ----------
    points : array-like
        The nxm array of points in the line, with m >= 2.
    crossings : array-like or None
        The nx2 array of crossing positions. If None, crossings
        are not plotted. Defaults to None.
    '''
    import matplotlib.pyplot as plt

    if fig_ax is not None:
        fig, ax = fig_ax
    else:
        fig, ax = plt.subplots()
    ax.plot(points[:, 0], points[:, 1])
    ax.set_xticks([])
    ax.set_yticks([])

    xmin, ymin = n.min(points[:, :2], axis=0)
    xmax, ymax = n.max(points[:, :2], axis=0)
    dx = (xmax - xmin) / 10.
    dy = (ymax - ymin) / 10.

    ax.set_xlim(xmin - dx, xmax + dx)
    ax.set_ylim(ymin - dy, ymax + dy)

    if mark_start:
        ax.plot([points[0, 0]], [points[0, 1]], color='blue',
                marker='o')

    if crossings is not None and len(crossings):
        crossings = n.array(crossings)
        ax.plot(crossings[:, 0], crossings[:, 1], 'ro', alpha=0.5)
    if show:
        fig.show()

    return fig, ax

def plot_shell(func, points, mode='auto', **kwargs):
    if mode =='auto':
        try:
            import vispy
            mode = 'vispy'
        except ImportError:
            pass
    if mode == 'auto':
        try:
            import mayavi.mlab as may
            mode = 'mayavi'
        except ImportError:
            pass

    if mode == 'mayavi':
        plot_shell_mayavi(func, points, **kwargs)
    elif mode == 'vispy':
        plot_shell_vispy(func, points, **kwargs)
    else:
        raise ValueError('invalid toolkit/mode')

def plot_shell_mayavi(func,
                      points,
                      number_of_samples=10,
                      zero_centroid=False,
                      sphere_radius_factor=2.,
                      opacity=0.3, **kwargs):
    positions, values = func(
        number_of_samples, zero_centroid=zero_centroid)

    thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi / 2.
    phis = n.linspace(0, 2 * n.pi, 157)

    thetas, phis = n.meshgrid(thetas, phis)

    r = sphere_radius_factor * n.max(points)
    zs = r * n.cos(thetas)
    xs = r * n.sin(thetas) * n.cos(phis)
    ys = r * n.sin(thetas) * n.sin(phis)

    import mayavi.mlab as may

    may.mesh(xs, ys, zs, scalars=values, opacity=opacity, **kwargs)

def plot_sphere_shell_vispy(func, rows=100, cols=100,
                            radius=1.,
                            opacity=0.5,
                            method='latitude',
                            edge_color=None,
                            cmap='hsv',
                            smooth=0,
                            **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    from vispy.scene import Sphere, ArcballCamera

    s = Sphere(rows=rows, cols=cols, method=method,
               edge_color=edge_color,
               radius=radius)
    mesh = s.mesh
    md = mesh._meshdata
    vertices = md.get_vertices()

    values = n.zeros(len(vertices))

    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)
        vertex = vertex / n.sqrt(n.sum(vertex*vertex))
        theta = n.arccos(vertex[2])
        phi = n.arctan2(vertex[1], vertex[0])

        if n.isnan(theta):
            theta = 0.0
        values[i] = func(theta=theta, phi=phi)
    vprint()

    colours = n.zeros((len(values), 4))
    max_val = n.max(values)
    min_val = n.min(values)
    unique_values = n.unique(colours)
    max_val += (1. + 1./len(unique_values))*(max_val - min_val)
    diff = (max_val - min_val)

    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    for i in range(len(colours)):
        colours[i] = cm(((values[i] - min_val) / diff))

    colours[:, -1] = opacity

    faces = md.get_faces()
    for si in range(smooth):
        new_colours = [[n.array(row) for _ in range(3)]
                       for row in colours]
        for i, face in enumerate(faces):
            new_colours[face[0]].append(colours[face[1]])
            new_colours[face[0]].append(colours[face[2]])
            new_colours[face[1]].append(colours[face[0]])
            new_colours[face[1]].append(colours[face[2]])
            new_colours[face[2]].append(colours[face[0]])
            new_colours[face[2]].append(colours[face[1]])

        new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])

        colours = new_colours

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    vispy_canvas.view.camera = ArcballCamera(fov=30, distance=radius*6)
    vispy_canvas.view.add(s)
    vispy_canvas.show()

def plot_colours_sphere_shell_vispy(func, rows=100, cols=100,
                                    radius=1.,
                                    opacity=0.6,
                                    method='latitude',
                                    edge_color=None,
                                    smooth=0,
                                    virtual=True,
                                    jones=True,
                                    **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    from vispy.scene import Sphere, ArcballCamera
    import csv

    s = Sphere(rows=rows, cols=cols, method=method,
               edge_color=edge_color,
               radius=radius)
    mesh = s.mesh
    md = mesh._meshdata
    vertices = md.get_vertices()

    values = []
    #results_file = open('colours_and_angles_sphere' + str(virtual) + '.csv', 'w')
    #results_file.close()

    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)
        vertex = vertex / n.sqrt(n.sum(vertex*vertex))
        theta = n.arccos(vertex[2])
        phi = n.arctan2(vertex[1], vertex[0])

        if n.isnan(theta):
            theta = 0.0

        try:
            with Timeout(300):
                value = func(theta=theta, phi=phi, virtual=virtual, jones=jones)
                values.append(value)
        except Timeout.Timeout:
            values.append(values[-1])

    vprint()

    colours = values

    colours = [[x[0], x[1], x[2], opacity] for x in colours]

    faces = md.get_faces()
    for si in range(smooth):
        new_colours = [[n.array(row) for _ in range(3)]
                       for row in colours]
        for i, face in enumerate(faces):
            new_colours[face[0]].append(colours[face[1]])
            new_colours[face[0]].append(colours[face[2]])
            new_colours[face[1]].append(colours[face[0]])
            new_colours[face[1]].append(colours[face[2]])
            new_colours[face[2]].append(colours[face[0]])
            new_colours[face[2]].append(colours[face[1]])

        new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])

        colours = new_colours

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    vispy_canvas.view.camera = ArcballCamera(fov=30, distance=radius*6)
    vispy_canvas.view.add(s)
    vispy_canvas.show()


def plot_shell_vispy(func,
                     points,
                     number_of_samples=10,
                     radius=None,
                     zero_centroid=False,
                     sphere_radius_factor=2.,
                     opacity=0.5,
                     cmap='hsv',
                     **kwargs):
    '''func must be a function returning values at angles and points,
    like OpenKnot._alexander_map_values.
    '''

    positions, values = func(
        number_of_samples, radius=radius,
        zero_centroid=zero_centroid)

    thetas = n.arcsin(n.linspace(-1, 1, 100)) + n.pi/2.
    phis = n.linspace(0, 2*n.pi, 157)

    thetas, phis = n.meshgrid(thetas, phis)

    r = sphere_radius_factor*n.max(points)
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
    cm = plt.get_cmap(cmap)
    for i in range(colours.shape[0]):
        for j in range(colours.shape[1]):
            colours[i, j] = cm(((values[i, j] - min_val) / diff))

    colours[:, :, -1] = opacity

    from vispy.scene import GridMesh
    from pyknot2.visualise import vispy_canvas
    mesh = GridMesh(xs, ys, zs, colors=colours)
    vispy_canvas.view.add(mesh)


def plot_cell(lines, mode='auto', **kwargs):
    if mode == 'auto':
        try:
            import vispy
            mode = 'vispy'
        except ImportError:
            pass
    if mode == 'auto':
        try:
            import mayavi.mlab as may
            mode = 'mayavi'
        except ImportError:
            pass

    if mode == 'mayavi':
        plot_cell_mayavi(lines, **kwargs)
    elif mode == 'vispy':
        plot_cell_vispy(lines, zero_centroid=False, **kwargs)
    else:
        raise ValueError('invalid toolkit/mode')

def plot_cell_mayavi(lines, boundary=None, clf=True, smooth=True,
                     min_seg_length=5, **kwargs):
    import mayavi.mlab as may
    may.clf()

    hues = n.linspace(0, 1, len(lines) + 1)[:-1]
    colours = [hsv_to_rgb(hue, 1, 1) for hue in hues]
    random.shuffle(colours)
    i = 0
    for (line, colour) in zip(lines, colours):
        vprint('Plotting line {} / {}\r'.format(i, len(lines)-1),
               False)
        i += 1
        for segment in line:
            if len(segment) < min_seg_length:
                continue
            plot_line(segment, mode='mayavi',
                      clf=False, color=colour, **kwargs)

    if boundary is not None:
        draw_bounding_box_mayavi(boundary, **kwargs)

def draw_bounding_box_mayavi(shape, colour=(0, 0, 0), tube_radius=1, markz=False):
    if shape is not None:
        if isinstance(shape, (float, int)):
            shape = ensure_shape_tuple(shape)
        if len(shape) == 3:
            shape = (0, shape[0], 0, shape[1], 0, shape[2])
    import mayavi.mlab as may

    xmin, xmax, ymin, ymax, zmin, zmax = shape
    ls = []
    ls.append(n.array([[xmax, ymax, zmin],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmax, ymin, zmin],[xmax, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmin],[xmin, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmin, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmax],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmax],[xmax, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmin],[xmax, ymax, zmin]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmax, ymin, zmin]]))
    ls.append(n.array([[xmax, ymin, zmax],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmax],[xmin, ymax, zmax]]))
    ls.append(n.array([[xmax, ymin, zmin],[xmax, ymax, zmin]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmin, ymax, zmin]]))

    ls = [interpolate(p) for p in ls]

    for line in ls:
        may.plot3d(line[:, 0], line[:, 1], line[:, 2],
                   color=colour, tube_radius=tube_radius)

def plot_cell_vispy(lines, boundary=None, clf=True, colours=None,
                    randomise_colours=True, tube_radius=1., **kwargs):
    if clf:
        clear_vispy_canvas()

    if colours is None:
        hues = n.linspace(0, 1, len(lines) + 1)[:-1]
        colours = [hsv_to_rgb(hue, 1, 1) for hue in hues]
        if randomise_colours:
            random.shuffle(colours)
    i = 0
    segments = []
    segment_colours = []
    segment_radii = []
    if not isinstance(tube_radius, list):
        tube_radius = [tube_radius for _ in range(len(lines))]
    for (line, colour, radius) in zip(lines, colours, tube_radius):
        vprint('Plotting line {} / {}\r'.format(i, len(lines)-1),
               False)
        i += 1
        for segment in line:
            if len(segment) < 4:
                continue
            segments.append(segment)
            segment_colours.append(colour)
            segment_radii.append(radius)
            # plot_line_vispy(segment,
            #                 clf=False, colour=colour, **kwargs)
    plot_lines_vispy(segments, colours=segment_colours,
                     tube_radius=segment_radii, **kwargs)

    if boundary is not None:
        draw_bounding_box_vispy(boundary, tube_radius=tube_radius[0])


def draw_bounding_box_vispy(shape, colour=(0, 0, 0), tube_radius=1.):
    if shape is not None:
        if isinstance(shape, (float, int)):
            shape = ensure_shape_tuple(shape)
        if len(shape) == 3:
            shape = (0, shape[0], 0, shape[1], 0, shape[2])

    xmin, xmax, ymin, ymax, zmin, zmax = shape
    ls = []
    ls.append(n.array([[xmax, ymax, zmin],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmax, ymin, zmin],[xmax, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmin],[xmin, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmin, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmax],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmax],[xmax, ymin, zmax]]))
    ls.append(n.array([[xmin, ymax, zmin],[xmax, ymax, zmin]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmax, ymin, zmin]]))
    ls.append(n.array([[xmax, ymin, zmax],[xmax, ymax, zmax]]))
    ls.append(n.array([[xmin, ymin, zmax],[xmin, ymax, zmax]]))
    ls.append(n.array([[xmax, ymin, zmin],[xmax, ymax, zmin]]))
    ls.append(n.array([[xmin, ymin, zmin],[xmin, ymax, zmin]]))

    ls = [interpolate(p) for p in ls]

    # plot_lines_vispy(ls, colours=['black' for _ in ls],
    #                  tube_radius=tube_radius, zero_centroid=False)
    plot_lines_vispy(ls, clf=False, colours=[colour for _ in ls],
                     zero_centroid=False, tube_radius=tube_radius)
    # for line in ls:
    #     plot_line_vispy(line, clf=False, colour=colour,
    #                     zero_centroid=False, tube_radius=tube_radius)

    global vispy_canvas
    vispy_canvas.central_widget.children[0].camera.center = (
        shape[1] / 2., shape[1] / 2., shape[1] / 2.)
    vispy_canvas.update()

def interpolate(p, num=10):
    p1, p2 = p
    return n.array([n.linspace(i, j, num) for i, j in zip(p1, p2)]).T

def cell_to_povray(filen, lines, shape):
    from jinja2 import Environment, FileSystemLoader
    env = Environment(loader=FileSystemLoader(
        '/home/asandy/devel/pyknot2/pyknot2/templates'))
    template = env.get_template('cell.pov')

    colours = n.linspace(0, 1, len(lines) + 1)[:-1]
    colours = n.array([hsv_to_rgb(c, 1, 1) for c in colours])

    coloured_segments = []
    for line, colour in zip(lines, colours):
        for segment in line:
            if len(segment) > 3:
                coloured_segments.append((segment, colour))

    with open(filen, 'w') as fileh:
        fileh.write(template.render(lines=coloured_segments))


def vispy_save_png(filename, region=None, size=None):
    img = vispy_canvas.render(region, size)
    import vispy.io as io
    io.write_png(filename, img)

def plot_sphere_mollweide_vispy(func, circle_points=50, depth=2,
                                edge_color=None,
                                cmap='hsv',
                                smooth=0,
                                mesh='circles',
                                **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    # from vispy.scene import Sphere, ArcballCamera
    from vispy.scene import TurntableCamera, Mesh

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)
    vertices[:, 0] *= 2*n.sqrt(2)
    vertices[:, 1] *= n.sqrt(2)
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    values = n.zeros(len(vertices))

    thetas = []
    phis = []
    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)

        intermediate = n.arcsin(vertex[1] / n.sqrt(2))
        theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)
        phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))

        # theta = n.arccos(vertex[2])
        # phi = n.arctan2(vertex[1], vertex[0])

        if n.isnan(theta):
            theta = 0.0
            print('theta', vertex)
        if n.isnan(phi):
            phi = 0.0
            print('phi', vertex)

        thetas.append(theta)
        phis.append(phi)
        values[i] = func(theta + n.pi/2, phi + n.pi)
    vprint()

    print('thetas', n.min(thetas), n.max(thetas))
    print('phis', n.min(phis), n.max(phis))

    colours = n.zeros((len(values), 4))
    max_val = n.max(values)
    min_val = n.min(values)
    unique_values = n.unique(colours)
    max_val += (1. + 1./len(unique_values))*(max_val - min_val)
    diff = (max_val - min_val)

    import matplotlib.pyplot as plt
    cm = plt.get_cmap(cmap)
    for i in range(len(colours)):
        colours[i] = cm(((values[i] - min_val) / diff))

    colours[:, -1] = 1.

    faces = md.get_faces()
    for si in range(smooth):
        new_colours = [[n.array(row) for _ in range(3)]
                       for row in colours]
        for i, face in enumerate(faces):
            new_colours[face[0]].append(colours[face[1]])
            new_colours[face[0]].append(colours[face[2]])
            new_colours[face[1]].append(colours[face[0]])
            new_colours[face[1]].append(colours[face[2]])
            new_colours[face[2]].append(colours[face[0]])
            new_colours[face[2]].append(colours[face[1]])

        new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])

        colours = new_colours

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    clear_vispy_canvas()
    vispy_canvas.view.camera = TurntableCamera(fov=1, distance=300, elevation=90., azimuth=0.)
    vispy_canvas.view.add(mesh)
    vispy_canvas.show()

def plot_colours_sphere_mollweide_vispy(func, circle_points=50, depth=2,
                                        smooth=0,
                                        mesh='circles',
                                        virtual=True,
                                        jones=False,
                                        mode='lambert',
                                        **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    # from vispy.scene import Sphere, ArcballCamera
    from vispy.scene import TurntableCamera, Mesh
    import csv

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)

    if mode == 'lambert':
        vertices[:, 0] *= 2
        vertices[:, 1] *= 2
    else:
        vertices[:, 0] *= 2*n.sqrt(2)
        vertices[:, 1] *= n.sqrt(2)
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    values = []
    #results_file = open('colours_and_angles_map' + str(virtual) + '.csv', 'w')
    #results_file.close()

    thetas = []
    phis = []
    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)

        if mode == 'lambert':
            rho = n.sqrt(vertex[0]**2 + vertex[1]**2)
            theta = 2 * n.arccos(rho / 2)
            phi = n.arctan2(vertex[1], vertex[0])

        else:
            intermediate = n.arcsin(vertex[1] / n.sqrt(2))
            theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)
            phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))

        if n.isnan(theta):
            theta = 0.0
            print('theta', vertex)
        if n.isnan(phi):
            phi = 0.0
            print('phi', vertex)

        thetas.append(theta)
        phis.append(phi)
        if mode == 'lambert':
            try:
                with Timeout(300):
                    value = func(theta=theta, phi=phi, virtual=virtual, jones=jones)
                    values.append(value)
            except Timeout.Timeout:
                values.append(values[-1])

        else:
            value = func(theta=theta + n.pi/2, phi=phi + n.pi, virtual=virtual, jones=jones)
            values.append(value)
        #with open('colours_and_angles_map' + str(virtual) + '.csv', 'a') as results_file:
        #    writer = csv.writer(results_file)
        #    writer.writerow(value)

    vprint()

    return thetas, phis, values

    print('thetas', n.min(thetas), n.max(thetas))
    print('phis', n.min(phis), n.max(phis))

    colours = values

    colours = [[x[0], x[1], x[2], 1.] for x in colours]

    faces = md.get_faces()
    for si in range(smooth):
        new_colours = [[n.array(row) for _ in range(3)]
                       for row in colours]
        for i, face in enumerate(faces):
            new_colours[face[0]].append(colours[face[1]])
            new_colours[face[0]].append(colours[face[2]])
            new_colours[face[1]].append(colours[face[0]])
            new_colours[face[1]].append(colours[face[2]])
            new_colours[face[2]].append(colours[face[0]])
            new_colours[face[2]].append(colours[face[1]])

        new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])

        colours = new_colours

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    clear_vispy_canvas()
    vispy_canvas.view.camera = TurntableCamera(fov=1, distance=350, elevation=90., azimuth=0.)
    vispy_canvas.view.add(mesh)
    vispy_canvas.show()

def plot_sphere_lambert_sharp_vispy(func, circle_points=50, depth=2, virtual=True, jones=True,
                                    output_size=500, mesh='circles', mode='mollweide',
                                    **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    from vispy.scene import TurntableCamera, Mesh

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)
    if mode == 'lambert':
        vertices[:, 0] *= 2
        vertices[:, 1] *= 2
    else:
        vertices[:, 0] *= 2*n.sqrt(2)
        vertices[:, 1] *= n.sqrt(2)
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    values = []

    thetas = []
    phis = []
    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)

        if mode == 'lambert':
            theta = 2*n.arccos(n.sqrt(vertex[0]**2 + vertex[1]**2)/2.)
            phi = n.arctan2(vertex[1], vertex[0])
        else:
            intermediate = n.arcsin(vertex[1] / n.sqrt(2))
            theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)
            phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))

        if n.isnan(theta):
            theta = 0.0
            print('theta', vertex)
        if n.isnan(phi):
            phi = 0.0
            print('phi', vertex)

        thetas.append(theta)
        phis.append(phi)
        if mode == 'lambert':
            try:
                with Timeout(300):
                    value = func(theta=theta, phi=phi, virtual=virtual, jones=jones)
                    values.append(value)
            except Timeout.Timeout:
                values.append(values[-1])

        else:
            try:
                with Timeout(300):
                    value = func(theta=theta + n.pi/2, phi=phi + n.pi, virtual=virtual, jones=jones)
                    values.append(value)
            except Timeout.Timeout:
                values.append(values[-1])
        #
        # try:
        #     with Timeout(300):
        #         value = func(theta=theta, phi=phi, virtual=virtual, jones=jones)
        #         values.append(value)
        # except Timeout.Timeout:
        #     values.append(values[-1])
    vprint()

    import svgwrite as svg
    d = svg.Drawing()

    for tri_i, triangle in enumerate(indices):
        v1 = vertices[triangle[0]]
        v2 = vertices[triangle[1]]
        v3 = vertices[triangle[2]]

        c1 = [int(x*255) for x in list(values[triangle[0]])]
        c2 = [int(x*255) for x in list(values[triangle[1]])]
        c3 = [int(x*255) for x in list(values[triangle[2]])]

        if mode == 'lambert':
            offset_x = 2.001412
        else:
            offset_x = 2.001412 * n.sqrt(2)
        offset_y = 2.0214
        v1_x = (v1[0] + offset_x) / 4. * output_size
        v1_y = (v1[1] + offset_y) / 4. * output_size

        v2_x = (v2[0] + offset_x) / 4. * output_size
        v2_y = (v2[1] + offset_y) / 4. * output_size

        v3_x = (v3[0] + offset_x) / 4. * output_size
        v3_y = (v3[1] + offset_y) / 4. * output_size

        c_x = n.average([v1_x, v2_x, v3_x])
        c_y = n.average([v1_y, v2_y, v3_y])

        v12_x = (v1_x + v2_x) / 2.
        v12_y = (v1_y + v2_y) / 2.
        v13_x = (v1_x + v3_x) / 2.
        v13_y = (v1_y + v3_y) / 2.
        v23_x = (v2_x + v3_x) / 2.
        v23_y = (v2_y + v3_y) / 2.

        d.add(d.polygon([[v1_x, v1_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v13_x, v13_y]],
                        fill='rgb({},{},{})'.format(*c1[:3]),
                        stroke='rgb({},{},{})'.format(*c1[:3]),
                        stroke_width='0.5'))

        d.add(d.polygon([[v2_x, v2_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*c2[:3]),
                        stroke='rgb({},{},{})'.format(*c2[:3]),
                        stroke_width='0.5'))

        d.add(d.polygon([[v3_x, v3_y],
                         [v13_x, v13_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*c3[:3]),
                        stroke='rgb({},{},{})'.format(*c3[:3]),
                        stroke_width='0.5'))

    return d

# run 'd.saveas('filename.svg') to save the result. You can view the svg using your browser.


def circles_ellipse_mesh(radial=5, azimuthal=100):
    angles = n.linspace(0, 2*n.pi, azimuthal + 1)[:-1]
    radii = n.linspace(0, 1., radial)

    offsets = n.zeros(len(angles))
    offsets[::2] += 2*n.pi / azimuthal

    vertices = []
    for i, radius in enumerate(radii):
        cur_angles = angles
        if i % 2 == 0:
            cur_angles += 0.5 * (2*n.pi) / azimuthal
        points = n.zeros((len(angles), 3))
        points[:, 0] = n.cos(angles) * radius
        points[:, 1] = n.sin(angles) * radius

        vertices.append(points)

    vertices = n.vstack(vertices)

    indices = []
    num_angles = len(angles)
    for num, radius in enumerate(radii[:-1]):
        base_index = num_angles * num
        next_num = num + 1
        for i in range(len(angles)):
            cur_index = base_index + i
            next_index = base_index + (i + 1) % num_angles
            next_r_index_1 = base_index + ((i - 1) % num_angles) + num_angles
            next_r_index_2 = base_index + i + num_angles
            indices.append((cur_index, next_r_index_1, next_r_index_2))
            indices.append((cur_index, next_r_index_2, next_index))

    return vertices, n.array(indices)

def plot_vertices_indices(vertices, indices):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    ax.plot(vertices[:, 0], vertices[:, 1], 'o')

    for triangle in indices:
        ax.plot([vertices[triangle[0]][0],
                 vertices[triangle[1]][0],
                 vertices[triangle[2]][0]],
                [vertices[triangle[0]][1],
                 vertices[triangle[1]][1],
                 vertices[triangle[2]][1]])

    fig.show()
    return fig, ax

def recursive_ellipse_mesh(init_num_points=10, depth=2):
    triangle = n.array([[0., 0.],
                        [1., 0.],
                        [0.5, 1.5]])

    angles = n.linspace(0, 2*n.pi, init_num_points + 1)[:-1]
    xs = n.sin(angles)
    ys = n.cos(angles)
    centre = n.array([0., 0.])

    vertices = []
    indices = []
    vertices_dict = {}

    for i in range(len(angles)):
        triangle = n.array([centre, (xs[i], ys[i]),
                            (xs[(i+1) % len(angles)],
                             ys[(i+1) % len(angles)])])
        indices.extend(get_subtriangles(triangle, vertices,
                                        vertices_dict, 0,
                                        depth))
    # indices = get_subtriangles(triangle, vertices, {}, 0, depth)

    return (n.array(vertices), n.array(indices))

def get_subtriangles(points, vertices, vertices_dict, depth, max_depth):
    if depth > max_depth:
        indices = []
        for point in points:
            point = tuple(point)
            if tuple(point) not in vertices_dict:
                vertices.append(point)
                vertices_dict[point] = len(vertices) - 1
            indices.append(vertices_dict[point])

        return [indices]

    arr_points = n.array(points)

    centre = n.average(arr_points, axis=0)
    half_edges = arr_points + 0.5*(n.roll(arr_points, -1, axis=0) - arr_points)

    new_points = n.array([arr_points[0],
                          half_edges[0],
                          arr_points[1],
                          half_edges[1],
                          arr_points[2],
                          half_edges[2]])

    subtriangles = []
    for i in range(6):
        new_triangle = n.array([centre, new_points[i], new_points[(i+1) % 6]])
        subtriangles.extend(get_subtriangles(new_triangle, vertices,
                                             vertices_dict,
                                             depth+1, max_depth))
    return subtriangles

def canvas_size(x=800, y=600):
    ensure_vispy_canvas()
    canvas = vispy_canvas
    from vispy import app, scene, color

    canvas.view.resize(x, y)
    canvas.show()

import numpy as np

def plot_sphere_lambert_sharp_vispy_sandy(func, circle_points=50, depth=2,
                                    output_size=500,
                                    edge_color=None,
                                    cmap='brg',
                                    smooth=0,
                                    mesh='circles',
                                    **kwargs):
    '''func must be a function of sphere angles theta, phi'''

    from vispy.scene import TurntableCamera, Mesh

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)
    vertices[:, 0] *= 2*n.sqrt(2)  # changed for mollweide
    vertices[:, 1] *= n.sqrt(2)  # changed for mollweide
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    values = n.zeros(len(vertices))

    print('pre')
    thetas = []
    phis = []
    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)

        intermediate = n.arcsin(vertex[1] / n.sqrt(2))  # mollweide
        theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)  # mollweide
        phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))  # mollweide

        if n.isnan(theta):
            theta = 0.0
            print('theta', vertex)
        if n.isnan(phi):
            phi = 0.0
            print('phi', vertex)

        thetas.append(theta)
        phis.append(phi)
        values[i] = func(theta + n.pi/2, phi + n.pi)
    vprint()

    import svgwrite as svg
    d = svg.Drawing()

    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(cmap)

    max_value = n.max(values)
    min_value = n.min(values)
    def normalise(v):
        return (v - min_value) / (max_value - min_value)
    for tri_i, triangle in enumerate(indices):
        v1 = vertices[triangle[0]]
        v2 = vertices[triangle[1]]
        v3 = vertices[triangle[2]]

        c1 = values[triangle[0]]
        c2 = values[triangle[1]]
        c3 = values[triangle[2]]

        offset_x = 2.001412 * 1.5 #changed for mollweide
        offset_y = 2.0214
        v1_x = (v1[0] + offset_x) / 4. * output_size
        v1_y = (v1[1] + offset_y) / 4. * output_size

        v2_x = (v2[0] + offset_x) / 4. * output_size
        v2_y = (v2[1] + offset_y) / 4. * output_size

        v3_x = (v3[0] + offset_x) / 4. * output_size
        v3_y = (v3[1] + offset_y) / 4. * output_size

        c_x = n.average([v1_x, v2_x, v3_x])
        c_y = n.average([v1_y, v2_y, v3_y])

        v12_x = (v1_x + v2_x) / 2.
        v12_y = (v1_y + v2_y) / 2.
        v13_x = (v1_x + v3_x) / 2.
        v13_y = (v1_y + v3_y) / 2.
        v23_x = (v2_x + v3_x) / 2.
        v23_y = (v2_y + v3_y) / 2.

        cmap_1 = cmap(normalise(c1))
        rgb_cmap1 = [int(c*255) for c in cmap_1[:3]]
        d.add(d.polygon([[v1_x, v1_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v13_x, v13_y]],
                        fill='rgb({},{},{})'.format(*rgb_cmap1),
                        stroke='rgb({},{},{})'.format(*rgb_cmap1),
                        stroke_width='0.5'))

        cmap_2 = cmap(normalise(c2))
        rgb_cmap2 = [int(c*255) for c in cmap_2[:3]]
        d.add(d.polygon([[v2_x, v2_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*rgb_cmap2),
                        stroke='rgb({},{},{})'.format(*rgb_cmap2),
                        stroke_width='0.5'))

        cmap_3 = cmap(normalise(c3))
        rgb_cmap3 = [int(c*255) for c in cmap_3[:3]]
        d.add(d.polygon([[v3_x, v3_y],
                         [v13_x, v13_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*rgb_cmap3),
                        stroke='rgb({},{},{})'.format(*rgb_cmap3),
                        stroke_width='0.5'))

        # print('cols', c1, cmap_1, c2, cmap_2, c3, cmap_3)

    return d

def mollweide_values(curve, circle_points=159, depth=316, mode='mollweide',
                     mesh='circles', virtual=True, jones=True):
    '''func must be a function of sphere angles theta, phi'''

    from vispy.scene import Mesh

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)
    if mode == 'lambert':
        vertices[:, 0] *= 2
        vertices[:, 1] *= 2
    else:
        vertices[:, 0] *= 2*n.sqrt(2)
        vertices[:, 1] *= n.sqrt(2)
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    values = []

    thetas = []
    phis = []
    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)

        if mode == 'lambert':
            theta = 2*n.arccos(n.sqrt(vertex[0]**2 + vertex[1]**2)/2.)
            phi = n.arctan2(vertex[1], vertex[0])
        else:
            intermediate = n.arcsin(vertex[1] / n.sqrt(2))
            theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)
            phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))

        if n.isnan(theta):
            theta = 0.0
            print('theta', vertex)
        if n.isnan(phi):
            phi = 0.0
            print('phi', vertex)

        thetas.append(theta)
        phis.append(phi)
        if mode == 'lambert':
            try:
                with Timeout(300):
                    value = curve.invariants_to_colours(theta=theta, phi=phi, virtual=virtual, jones=jones)
                    values.append(value)
            except Timeout.Timeout:
                values.append(values[-1])
        else:
            try:
                with Timeout(300):
                    value = curve.invariants_to_colours(theta=theta + n.pi/2, phi=phi + n.pi, virtual=virtual, jones=jones)
                    values.append(value)
            except Timeout.Timeout:
                values.append(values[-1])
        # try:
        #     with Timeout(300):
        #         value = curve.invariants_to_colours(theta=theta, phi=phi, virtual=virtual, jones=jones)
        #         values.append(value)
        # except Timeout.Timeout:
        #     values.append(values[-1])

    return thetas, phis, values

def sphere_values(curve, rows=159, cols=316, method='latitude', virtual=True, jones=True):

    from vispy.geometry import create_sphere

    mesh = create_sphere(rows=rows, cols=cols)

    vertices = mesh.get_vertices()

    values = []
    thetas = []
    phis = []

    for i, vertex in enumerate(vertices):
        if i % 10 == 0:
            vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)
        vertex = vertex / n.sqrt(n.sum(vertex*vertex))
        theta = n.arccos(vertex[2])
        phi = n.arctan2(vertex[1], vertex[0])

        if n.isnan(theta):
            theta = 0.0

        thetas.append(theta)
        phis.append(phi)

        try:
            with Timeout(300):
                value = curve.invariants_to_colours(theta=theta, phi=phi, virtual=virtual, jones=jones)
                values.append(value)
        except Timeout.Timeout:
            values.append(values[-1])

    return thetas, phis, values

def sphere_from_file(protein, source='/home/ka14927/Python_Files/map_values/', smooth=0,
                     virtual=True, low_res=False, radius=1, edge_colour=None):

    import json

    if virtual:
        inv_type = '_virtual_'
    else:
        inv_type = '_classical_'

    if low_res:
        filename = inv_type + 'sphere_values_10000'
        rows, cols = 71, 141
    else:
        filename = inv_type + 'sphere_values'
        rows, cols = 159, 316

    with open(source + protein + filename) as data_file:
        data = json.load(data_file)
    # thetas = [x[0] for x in data]
    # phis = [x[1] for x in data]
    values = [x[2] for x in data]

    from vispy.scene import Sphere, ArcballCamera

    # rows and cols mixed up between vispy.geometry.create_sphere and vispy.scene.Sphere
    s = Sphere(rows=cols, cols=rows, method='latitude',
               edge_color=edge_colour,
               radius=radius)
    mesh = s.mesh
    md = mesh._meshdata
    # vertices = md.get_vertices()

    colours = values

    colours = [[x[0], x[1], x[2], 1.0] for x in colours]

    faces = md.get_faces()
    for si in range(smooth):
        new_colours = [[n.array(row) for _ in range(3)]
                       for row in colours]
        for i, face in enumerate(faces):
            new_colours[face[0]].append(colours[face[1]])
            new_colours[face[0]].append(colours[face[2]])
            new_colours[face[1]].append(colours[face[0]])
            new_colours[face[1]].append(colours[face[2]])
            new_colours[face[2]].append(colours[face[0]])
            new_colours[face[2]].append(colours[face[1]])

        new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])

        colours = new_colours

    md.set_vertex_colors(colours)

    ensure_vispy_canvas()
    vispy_canvas.view.camera = ArcballCamera(fov=30, distance=radius*6)
    vispy_canvas.view.add(s)
    vispy_canvas.show()

def map_from_file(protein, source='/home/ka14927/Python_Files/map_values/', smooth=0,
                  virtual=True, low_res=False, mesh='circles', mode='mollweide', output_size=500):

    '''func must be a function of sphere angles theta, phi'''

    import json

    if virtual:
        inv_type = '_virtual_'
    else:
        inv_type = '_classical_'

    if low_res:
        filename = inv_type + 'map_values_10000'
        circle_points, depth = 71, 141
    else:
        filename = inv_type + 'map_values'
        circle_points, depth = 159, 316

    with open(source + protein + filename) as data_file:
        data = json.load(data_file)
    # thetas = [x[0] for x in data]
    # phis = [x[1] for x in data]
    values = [x[2] for x in data]

    from vispy.scene import TurntableCamera, Mesh

    if mesh == 'circles':
        vertices, indices = circles_ellipse_mesh(circle_points, depth)
    else:
        vertices, indices = recursive_ellipse_mesh(circle_points, depth)

    if mode == 'lambert':
        vertices[:, 0] *= 2
        vertices[:, 1] *= 2
    else:
        vertices[:, 0] *= 2*n.sqrt(2)
        vertices[:, 1] *= n.sqrt(2)
    mesh = Mesh(vertices, indices)

    md = mesh._meshdata
    vertices = md.get_vertices()

    colours = values

    colours = [[x[0], x[1], x[2], 1.] for x in colours]

    # # faces = md.get_faces()
    # # for si in range(smooth):
    # #     new_colours = [[n.array(row) for _ in range(3)]
    # #                    for row in colours]
    # #     for i, face in enumerate(faces):
    # #         new_colours[face[0]].append(colours[face[1]])
    # #         new_colours[face[0]].append(colours[face[2]])
    # #         new_colours[face[1]].append(colours[face[0]])
    # #         new_colours[face[1]].append(colours[face[2]])
    # #         new_colours[face[2]].append(colours[face[0]])
    # #         new_colours[face[2]].append(colours[face[1]])
    # #
    # #     new_colours = n.array([n.average(cs, axis=0) for cs in new_colours])
    # #
    # #     colours = new_colours
    # #
    # # md.set_vertex_colors(colours)
    # #
    # # ensure_vispy_canvas()
    # # clear_vispy_canvas()
    # # vispy_canvas.view.camera = TurntableCamera(fov=1, distance=350, elevation=90., azimuth=0.)
    # # vispy_canvas.view.add(mesh)
    # # vispy_canvas.show()
    # #
    # from vispy.scene import TurntableCamera, Mesh
    #
    # if mesh == 'circles':
    #     vertices, indices = circles_ellipse_mesh(circle_points, depth)
    # else:
    #     vertices, indices = recursive_ellipse_mesh(circle_points, depth)
    # if mode == 'lambert':
    #     vertices[:, 0] *= 2
    #     vertices[:, 1] *= 2
    # else:
    #     vertices[:, 0] *= 2*n.sqrt(2)
    #     vertices[:, 1] *= n.sqrt(2)
    # mesh = Mesh(vertices, indices)
    #
    # md = mesh._meshdata
    # vertices = md.get_vertices()
    #
    # values = []
    #
    # thetas = []
    # phis = []
    # for i, vertex in enumerate(vertices):
    #     if i % 10 == 0:
    #         vprint('\ri = {} / {}'.format(i, len(vertices)), newline=False)
    #
    #     if mode == 'lambert':
    #         theta = 2*n.arccos(n.sqrt(vertex[0]**2 + vertex[1]**2)/2.)
    #         phi = n.arctan2(vertex[1], vertex[0])
    #     else:
    #         intermediate = n.arcsin(vertex[1] / n.sqrt(2))
    #         theta = n.arcsin((2*intermediate + n.sin(2*intermediate)) / n.pi)
    #         phi = n.pi * vertex[0] / (2*n.sqrt(2) * n.cos(intermediate))
    #
    #
    #     if n.isnan(theta):
    #         theta = 0.0
    #         print('theta', vertex)
    #     if n.isnan(phi):
    #         phi = 0.0
    #         print('phi', vertex)
    #
    #     thetas.append(theta)
    #     phis.append(phi)
    #     try:
    #         with Timeout(300):
    #             value = func(theta=theta, phi=phi, virtual=virtual, jones=jones)
    #             values.append(value)
    #     except Timeout.Timeout:
    #         values.append(values[-1])
    # vprint()

    import svgwrite as svg
    d = svg.Drawing()

    for tri_i, triangle in enumerate(indices):
        v1 = vertices[triangle[0]]
        v2 = vertices[triangle[1]]
        v3 = vertices[triangle[2]]

        c1 = [int(x*255) for x in list(values[triangle[0]])]
        c2 = [int(x*255) for x in list(values[triangle[1]])]
        c3 = [int(x*255) for x in list(values[triangle[2]])]

        if mode == 'lambert':
            offset_x = 2.001412
        else:
            offset_x = 2.001412 * n.sqrt(2)
        offset_y = 2.0214
        v1_x = (v1[0] + offset_x) / 4. * output_size
        v1_y = (v1[1] + offset_y) / 4. * output_size

        v2_x = (v2[0] + offset_x) / 4. * output_size
        v2_y = (v2[1] + offset_y) / 4. * output_size

        v3_x = (v3[0] + offset_x) / 4. * output_size
        v3_y = (v3[1] + offset_y) / 4. * output_size

        c_x = n.average([v1_x, v2_x, v3_x])
        c_y = n.average([v1_y, v2_y, v3_y])

        v12_x = (v1_x + v2_x) / 2.
        v12_y = (v1_y + v2_y) / 2.
        v13_x = (v1_x + v3_x) / 2.
        v13_y = (v1_y + v3_y) / 2.
        v23_x = (v2_x + v3_x) / 2.
        v23_y = (v2_y + v3_y) / 2.

        d.add(d.polygon([[v1_x, v1_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v13_x, v13_y]],
                        fill='rgb({},{},{})'.format(*c1[:3]),
                        stroke='rgb({},{},{})'.format(*c1[:3]),
                        stroke_width='0.5'))

        d.add(d.polygon([[v2_x, v2_y],
                         [v12_x, v12_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*c2[:3]),
                        stroke='rgb({},{},{})'.format(*c2[:3]),
                        stroke_width='0.5'))

        d.add(d.polygon([[v3_x, v3_y],
                         [v13_x, v13_y],
                         [c_x, c_y],
                         [v23_x, v23_y]],
                        fill='rgb({},{},{})'.format(*c3[:3]),
                        stroke='rgb({},{},{})'.format(*c3[:3]),
                        stroke_width='0.5'))

    return d

