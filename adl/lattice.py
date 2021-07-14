import math


class square_lattice:
    def __init__(self, a, reps):
        self.a = a
        self.reps = reps
        self.coordinates = []
        for i in range(self.reps):
            for j in range(self.reps):
                x = i * self.a + self.a / 2
                y = j * self.a + self.a / 2
                self.coordinates.append((x, y))


class hexagon_lattice:
    def __init__(self, a, reps):
        self.a = a
        self.h = a * math.sin(60 * math.pi / 180)
        self.reps = reps
        self.coordinates = []
        # for i in range(self.reps):
        #     for j in range(self.reps + 1):
        #         x = i * self.a + self.a / 2
        #         y = j * self.a + self.a / 2
        #         if i % 2 == 1:
        #             y -= self.a / 2
        #         self.coordinates.append((x, y))


class honeycomb_lattice:
    def __init__(self, a, reps):
        self.a = a
        self.reps = reps
        self.coordinates = []
        for i in range(self.reps):
            for j in range(self.reps + 1):
                x = i * self.a + self.a / 2
                y = j * self.a + self.a / 2
                if i % 2 == 1:  # rows 1,3,5 ...
                    y -= self.a / 2
                    if j % 3 == 1:  # col 1,4,7 ...
                        continue
                if i % 2 == 0:  # rows 0,2,4 ...
                    if j % 3 == 2:  # col 2,5,8 ...
                        continue
                self.coordinates.append((x, y))
