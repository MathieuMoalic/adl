import math


class square_lattice:
    def __init__(self, adl):
        self.coordinates = [(adl.lattice_param / 2, adl.lattice_param / 2)]
        self.xsize = adl.lattice_param
        self.ysize = adl.lattice_param


class rectangular_lattice:
    def __init__(self, adl):
        self.coordinates = [(adl.lattice_param / 2, adl.lattice_param2 / 2)]
        self.xsize = adl.lattice_param
        self.ysize = adl.lattice_param2


class honeycomb_lattice:
    def __init__(self, adl):
        h = adl.lattice_param * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (adl.lattice_param, 0),
            (2 * adl.lattice_param, 0),
            (adl.lattice_param / 2, h),
            (2.5 * adl.lattice_param, h),
            (adl.lattice_param, 2 * h),
            (2 * adl.lattice_param, 2 * h),
        ]
        self.xsize = 3 * adl.lattice_param
        self.ysize = 2 * h


class hexagonal_lattice:
    def __init__(self, adl):
        h = adl.lattice_param * math.sin(60 * math.pi / 180)
        self.coordinates = [
            (0.25 * adl.lattice_param, 1.5 * h),
            (0.75 * adl.lattice_param, 0.5 * h),
            (
                1.25 * adl.lattice_param,
                1.5 * h,
            ),  # these extra 2 nodes prevent bad clipping
            (-0.25 * adl.lattice_param, 0.5 * h),  # needs top and bottom ones too
        ]
        self.xsize = adl.lattice_param
        self.ysize = 2 * h


class octagonal_lattice:
    def __init__(self, adl):
        self.coordinates = [
            (0.5 * adl.lattice_param, 0),
            (1.5 * adl.lattice_param, 0),
            (0.5 * adl.lattice_param, 2 * adl.lattice_param),
            (1.5 * adl.lattice_param, 2 * adl.lattice_param),
            (0, 0.5 * adl.lattice_param),
            (0, 1.5 * adl.lattice_param),
            (2 * adl.lattice_param, 0.5 * adl.lattice_param),
            (2 * adl.lattice_param, 1.5 * adl.lattice_param),
        ]
        self.xsize = 2 * adl.lattice_param
        self.ysize = 2 * adl.lattice_param
