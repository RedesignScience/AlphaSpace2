import nglview as nv


def view_structure(traj, surface_opacity=1.0, snapshot_idx=0):
    _ngl_view = nv.show_mdtraj(traj[snapshot_idx])

    if surface_opacity > 0:
        _ngl_view.add_surface(selection='protein', opacity=surface_opacity, color='white')

    return _ngl_view


def draw_sphere(view, coord, radius, color):
    view.shape.add_buffer("sphere", position=[coord], color=[color], radius=[radius])
    return True


def draw_cylinder(view, start, end, radius, color):
    view.shape.add_buffer("cylinder", position1=start, position2=end, color=color,
                          radius=radius)
    return True
