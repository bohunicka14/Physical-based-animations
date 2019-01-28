import tkinter
from tkinter import messagebox
from tkinter import *
import math


def multiply(m, array):
    """
    :param m:matrix
    :param array: list of coordinates
    :return: multiply list of points with matrix
    """
    for i in array:
        i = int(i)

    for i in range(0, len(array), 2):
        result_x = array[i] * m[0][0] + array[i + 1] * m[0][1]
        result_y = array[i] * m[1][0] + array[i + 1] * m[1][1]
        array[i] = result_x
        array[i + 1] = result_y
    return array


def unit_vector(vector):
    """
    :param vector
    :return: return vector with lenght = 1
    """
    length = math.sqrt(vector[0] ** 2 + vector[1] ** 2)
    return [vector[0] / length, vector[1] / length]


def projection(A, B, C):
    """
    :param A, B, C: vertex
    :return: projection of point C to line segment AB
    """
    x1 = A.x
    y1 = A.y
    x2 = B.x
    y2 = B.y
    x3 = C.x
    y3 = C.y

    px = x2 - x1
    py = y2 - y1
    dAB = px * px + py * py

    u = ((x3 - x1) * px + (y3 - y1) * py) / dAB

    x = x1 + u * px
    y = y1 + u * py
    return [x, y]


def mid_point(v1, v2):
    """
    :param v1, v2: vertex
    :return: middle point between two input vertices
    """
    v = [v2.x - v1.x, v2.y - v1.y]
    v[0], v[1] = v[0]/2, v[1]/2
    return [ v[0] + v1.x, v[1] + v1.y ]

class Feature():

    def mark(self):
        self.marked = True

    def clear(self):
        self.marked = False


class Vertex(Feature):

    def __init__(self, canvas, x, y, marked=False):
        """
        :param canvas: Canvas
        :param x: x point of vertex
        :param y: y point of vertex
        :param marked: True if its marked or False if its not
        """
        self.canvas = canvas
        self.x = x
        self.y = y
        self.marked = marked
        self.voronoi_region = [[self.x, self.y], [self.x, self.y]]
        self.vr1 = self.vr2 = None
        self.r = 5

    def update(self, x, y):
        """
        :param x: new x coordinate
        :param y: new y coordinate
        """
        self.x = x
        self.y = y
        self.voronoi_region = [[self.x, self.y], [self.x, self.y]]

    def draw_VR(self, color='blue', bold = False):
        """
        :param color: color of VR boundaries
        :param bold:
        :return: draws Voronoi plane
        """
        if self.vr2 is not None and self.vr1 is not None:
            self.canvas.delete(self.vr1)
            self.canvas.delete(self.vr2)
        if len(self.voronoi_region[0]) == 4:
            if bold:
                self.vr1 = self.canvas.create_line(self.voronoi_region[0], fill=color, width=3)
            else:
                self.vr1 = self.canvas.create_line(self.voronoi_region[0], fill=color)
        if len(self.voronoi_region[1]) == 4:
            if bold:
                self.vr2 = self.canvas.create_line(self.voronoi_region[1], fill=color, width=3)
            else:
                self.vr2 = self.canvas.create_line(self.voronoi_region[1], fill=color)

    def move_vr(self, x, y):
        """
        :param x: move VP for x points
        :param y: move VP for y points
        """
        self.canvas.move(self.vr1, x, y)
        self.canvas.move(self.vr2, x, y)

    def draw(self, canvas):
        """
        :param canvas: Canvas
        :return: draws vertex
        """
        self.canvas = canvas
        self.id = self.canvas.create_oval(self.x - self.r, self.y - self.r, self.x + self.r, self.y + self.r,
                                          fill='black')

    def update_color(self):
        """
        :return: draws marked vertex
        """
        self.id = self.canvas.create_oval(self.x - self.r, self.y - self.r, self.x + self.r, self.y + self.r,
                                          fill='red')

    def move(self, xx, yy):
        """
        :param xx: moves for xx points
        :param yy: moves for yy points
        :return: moves object
        """
        self.canvas.move(self.id, xx, yy)

    def __eq__(self, v):
        """
        :param v: vertex
        :return: return True if vertex is equivalent with input vertex
        """
        return self.x == v.x and self.y == v.y

    def __add__(self, v):
        """
        :param v: vertex
        :return: array, counts vertex with input vertex by their coordinates
        """
        return [self.x + v.x, self.y + v.y]

    def __sub__(self, v):
        """
        :param v: vertex
        :return: array, subtracts vertex with input vertex by their coordinates
        """
        return [self.x - v.x, self.y - v.y]

    def __mul__(self, k):
        """
        :param k: float
        :return: array, multiplies vertex by constant
        """
        return [self.x * k, self.y * k]


class Edge(Feature):

    def __init__(self, v1, v2, polygon, marked=False):
        """
        :param v1: first vertex
        :param v2: second vertex
        :param polygon: PolyObject reference
        :param marked: False if its not marked, True if its marked
        """
        self.v1 = v1
        self.v2 = v2
        self.polygon = polygon
        self.marked = marked
        self.voronoi_region = self.create_voronoi_region()
        self.vr1 = self.vr2 = None

    # def update_color(self):
    #     self.polygon.canvas.create_line(self.v1.x, self.v1.y, self.v2.x, self.v2.y, fill='orange',
    #                                     width=3)

    # def update(self, x1, y1, x2, y2):
    #     self.v1.update(x1, y1)
    #     self.v2.update(x2, y2)

    def get_directional_vector(self):
        """
        :return: directional vector
        """
        return [self.v2.x - self.v1.x, self.v2.y - self.v1.y]

    def get_normal_vector(self):
        """
        :return: mormal vector
        """
        directional = self.get_directional_vector()
        normal_vector = [-directional[1], directional[0]]
        vector = unit_vector(normal_vector)
        mid_point = [(self.v2.x - self.v1.x) / 2 + self.v1.x, (self.v2.y - self.v1.y) / 2 + self.v1.y]

        for edge in self.polygon.edges:
            if edge != self:
                # intersects = ray_line_segment_intersection(mid_point, vector, [edge.v1.x, edge.v1.y], [edge.v2.x, edge.v2.y])
                intersects = self.polygon.point_in_polygon(mid_point[0] + vector[0], mid_point[1] + vector[1])
                if intersects:
                    return [-vector[0], -vector[1]]
        return vector

    def get_inverse_normal_vector(self):
        """
        :return: inverse normal vector
        """
        normal = self.get_normal_vector()
        return [-normal[0], -normal[1]]

    def calculate_plain(self, v):
        """
        :param v: vertex
        :return: normal vector and constant of plain
        """
        n = unit_vector(self.get_directional_vector())

        if v == self.v2:
            n = [-n[0], -n[1]]

        c = n[0] * v.x + n[1] * v.y
        return [n, -c]

    def create_voronoi_region(self):
        """
        :return: properties of VR of edge
        """
        vector = self.get_normal_vector()

        vector[0] *= 1000
        vector[1] *= 1000

        ax = vector[0] + self.v1.x
        ay = vector[1] + self.v1.y

        bx = vector[0] + self.v2.x
        by = vector[1] + self.v2.y

        self.v1.voronoi_region[0].append(ax)
        self.v1.voronoi_region[0].append(ay)

        self.v2.voronoi_region[1].append(bx)
        self.v2.voronoi_region[1].append(by)

        return [[self.v1.x, self.v1.y, ax, ay], [self.v2.x, self.v2.y, bx, by]]

    def draw_VR(self, color='blue', bold = False):
        """
        :param color: color of VR boundaries
        :param bold:
        :return: draws Voronoi plane
        """
        if self.vr2 is not None and self.vr1 is not None:
            self.polygon.canvas.delete(self.vr1)
            self.polygon.canvas.delete(self.vr2)

        if bold:
            self.vr1 = self.polygon.canvas.create_line(self.voronoi_region[0], fill=color, width=3)
            self.vr2 = self.polygon.canvas.create_line(self.voronoi_region[1], fill=color, width=3)
        else:
            self.vr1 = self.polygon.canvas.create_line(self.voronoi_region[0], fill=color)
            self.vr2 = self.polygon.canvas.create_line(self.voronoi_region[1], fill=color)

    def move_vr(self, x, y):
        """
        :param x: moves for x points
        :param y: moves for y points
        """
        self.polygon.canvas.move(self.vr1, x, y)
        self.polygon.canvas.move(self.vr2, x, y)

    def update(self):
        self.voronoi_region = self.create_voronoi_region()

    def __eq__(self,e):
        """
        :param e: edge
        :return: True if edge is equivalent with input edge
        """
        return self.v1 == e.v1 and self.v2 == e.v2

    def __mul__(self, k):
        """
        :param k: float
        :return: multiplied directional vector of edge by constant
        """
        v = self.get_directional_vector()
        return [v[0]*k, v[1]*k]

class FeaturePair():

    def __init__(self, f1, f2):
        """
        :param f1: first feature, Edge/Vertex
        :param f2: second feature, Edge/Vertex
        """
        self.f1 = f1
        self.f2 = f2

    def type(self):
        """
        :return: type of case
        """
        result = ''
        if isinstance(self.f1, Edge):
            result += 'E'
        elif isinstance(self.f1, Vertex):
            result += 'V'

        if isinstance(self.f2, Edge):
            result += 'E'
        elif isinstance(self.f2, Vertex):
            result += 'V'

        return result

    def swap(self):
        """
        :return: swapped features f1, f2
        """
        tmp = self.f1
        self.f1 = self.f2
        self.f2 = tmp


class PolyObject:
    def __init__(self, canvas, x, y, tag):
        """
        :param canvas: Canvas
        :param x: x coordinate of polygon
        :param y: y coordinate of polygon
        :param tag: special tag for polygon
        """
        self.canvas = canvas
        self.x = x
        self.y = y
        self.tag = tag
        self.is_active = 0
        self.edges = []
        self.vertices = []
        self.voronoi_regions = []
        self.features = []

    def draw(self, coords, color):
        """
        :param coords: coordinates of polygon
        :param color: color of edges
        :return: draws a polygon
        """
        if self.is_active:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='red', width=3, tag=self.tag)
        else:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='black', width=1, tag=self.tag)

    def point_in_polygon(self, x, y):
        """
        :param x: integer
        :param y: integer
        :return: True if point lies in polygon
        """
        n = len(self.vertices)
        v = False

        b1x, b1y = self.vertices[0].x, self.vertices[0].y
        for i in range(n + 1):
            b2x, b2y = self.vertices[i % n].x, self.vertices[i % n].y
            if y > min(b1y, b2y):
                if y <= max(b1y, b2y):
                    if x <= max(b1x, b2x):
                        if b1y != b2y:
                            xin = (y - b1y) * (b2x - b1x) / (b2y - b1y) + b1x
                        if b1x == b2x or x <= xin:
                            v = not v
            b1x, b1y = b2x, b2y
        return v

    def set_active(self, state):
        """
        :param state: 1 or 0
        :return: changes outline of polygon
        """
        if state == 1:
            self.canvas.itemconfig(self.id, outline='red', width=3)
            self.is_active = 1
        else:
            self.canvas.itemconfig(self.id, outline='black', width=1)
            self.is_active = 0

    def draw_vertices(self):
        """
        :return: draws vertices of polygon
        """
        for vertex in self.vertices:
            vertex.draw(self.canvas)

    def set_features(self):
        """
        :return: array of vertices and edges of polygon
        """
        if len(self.vertices) == 0 and len(self.edges) == 0:
            ## init features
            for i in range(0, len(self.coords), 2):
                # print(len(self.coords), self.coords[i], self.coords[i+1])
                self.vertices.append(Vertex(self.canvas, self.coords[i], self.coords[i + 1]))
            self.draw_vertices()

            for i in range(len(self.vertices) - 1):
                self.edges.append(Edge(self.vertices[i], self.vertices[i + 1], self))
            self.edges.append(Edge(self.vertices[0], self.vertices[-1], self))
            # print(len(self.edges))
        else:
            ## update features according to new coords
            j = 0
            for i in range(0, len(self.coords), 2):
                self.vertices[j].update(self.coords[i], self.coords[i + 1])
                j += 1

            for e in self.edges:
                e.update()
        self.features = self.vertices + self.edges

    def set_coords(self, coords):
        """
        :param coords: coordinates of polygon
        """
        self.coords = coords
        self.set_features()

    def get_centroid(self):
        """
        :return: center of polygon
        """
        self.coords.append(self.coords[0])
        self.coords.append(self.coords[1])
        A = 0
        for i in range(0, len(self.coords) - 2, 2):
            A += self.coords[i] * self.coords[i + 3] - self.coords[i + 2] * self.coords[i + 1]
        A = A * (1 / 2)

        Cx = Cy = 0
        for i in range(0, len(self.coords) - 2, 2):
            Cx += (self.coords[i] + self.coords[i + 2]) * (
                    self.coords[i] * self.coords[i + 3] - self.coords[i + 2] * self.coords[i + 1])
            Cy += (self.coords[i + 1] + self.coords[i + 3]) * (
                    self.coords[i] * self.coords[i + 3] - self.coords[i + 2] * self.coords[i + 1])

        self.coords.pop()
        self.coords.pop()
        x = Cx / (6 * A)
        y = Cy / (6 * A)

        return [x, y]

    def rotate(self):
        """
        Calculates rotated coordinates of polygon
        """
        matrix = [[math.cos(math.radians(15)), math.sin(math.radians(15))],
                  [-math.sin(math.radians(15)), math.cos(math.radians(15))]]

        center = self.get_centroid()

        ## translate polygon to the center of coordinate system
        for i in range(0, len(self.coords), 2):
            self.coords[i] -= center[0]
            self.coords[i + 1] -= center[1]

        ## rotate polygon
        self.coords = multiply(matrix, self.coords)

        ## translate rotated polygon back to the original position
        for i in range(0, len(self.coords), 2):
            self.coords[i] += center[0]
            self.coords[i + 1] += center[1]

    def move(self, xx, yy):
        """
        :param xx: moves polygon for xx points
        :param yy: moves polygon for yy points
        """
        self.canvas.move(self.id, xx, yy)

    def __eq__(self,e):
        return self.v1 == e.v1 and self.v2 == e.v2


class Playground(Tk):
    def __init__(self):
        super().__init__()
        main = Frame(self, bg='red', width=1000)
        main.pack()

        self.id_segment = None
        self.highlighted_feature1 = None
        self.highlighted_feature2 = None

        self.labelframe = LabelFrame(main, text="Rotation", width=200)
        self.labelframe.pack(side=LEFT, fill=BOTH)

        rotate = Button(self.labelframe, text="Rotate", width=13, command=self.rotate_object)
        rotate.pack(side=TOP)

        self.playground = Canvas(main, width=800, height=600, highlightthickness=0)
        self.playground.pack(side=LEFT, fill=BOTH)
        self.update()

        self.polygons_array = []
        self.obj_id = None
        self.load_from_file('objects.txt')
        self.playground.bind("<ButtonPress-1>", self.click)
        self.playground.bind("<B1-Motion>", self.drag)
        self.playground.bind("<ButtonRelease-1>", self.drop)

    def rotate_object(self):
        """
        Rotates object
        """
        if self.p1.is_active:
            self.obj_id = self.p1.id
            self.p1.rotate()
            self.redraw_canvas(False)
            self.p1.set_coords(self.playground.coords(self.p1.id))
            for e in self.p1.edges:
                e.draw_VR('green')

            for e in self.p2.edges:
                e.draw_VR()


        elif self.p2.is_active:
            self.obj_id = self.p2.id
            self.p2.rotate()
            self.redraw_canvas(False)
            self.p2.set_coords(self.playground.coords(self.p2.id))
            for e in self.p1.edges:
                e.draw_VR('green')

            for e in self.p2.edges:
                e.draw_VR()

        else:
            messagebox.showwarning("Warning", "No object is selected!!!")
            return

    def load_from_file(self, file):
        """
        :param file: text file with coordinates
        """
        coords = []
        f = open(file, 'r')

        for line in f:
            coords.append([line.strip(' ')])
        f.close()
        self.create_polygons(coords)

    def create_polygons(self, polygons):
        """
        :param polygons: coordinates of polygons loaded from text file
        Creates and draws polygon
        """
        coords_1 = polygons[0]
        coords_2 = polygons[1]
        self.p1 = PolyObject(self.playground, 100, 100, 'object1')
        self.p2 = PolyObject(self.playground, 400, 400, 'object2')
        self.polygons_array.append(self.p1)
        self.polygons_array.append(self.p2)

        self.p1.draw(coords_1, 'white')
        self.p2.draw(coords_2, 'white')
        self.p1.set_coords(self.playground.coords(self.p1.id))
        self.p2.set_coords(self.playground.coords(self.p2.id))

    def redraw_canvas(self, draw_vr=True):
        """
        Redraws objects in canvas
        """
        self.playground.delete('all')
        self.p1.draw(self.p1.coords, 'white')
        self.p2.draw(self.p2.coords, 'white')

        if draw_vr:
            for e in self.p1.edges:
                e.draw_VR('green')

            for e in self.p2.edges:
                e.draw_VR()

    def click(self, event):
        """
        :param event: click event
        """
        self.obj_id = event.widget.find_closest(event.x, event.y)[0]
        if self.playground.gettags(event.widget.find_closest(event.x, event.y)):
            self.obj_tag = self.playground.gettags(event.widget.find_closest(event.x, event.y))[0]
        if 'object' in self.obj_tag:
            self.playground.lift(self.obj_id)
        if 'object1' in self.obj_tag:
            self.p1.set_active(1)
            self.p2.set_active(2)
        elif 'object2' in self.obj_tag:
            self.p2.set_active(1)
            self.p1.set_active(2)

        self.ex, self.ey = event.x, event.y

    def drop(self, event):
        """
        :param event: click event
        """
        if self.p1.is_active:
            self.p1.set_coords(self.playground.coords(self.p1.id))
        elif self.p2.is_active:
            self.p2.set_coords(self.playground.coords(self.p2.id))

    def drag(self, event):
        """
        :param event: click event
        """
        self.playground.move(self.obj_id, event.x - self.ex, event.y - self.ey)
        if self.p1.is_active:
            self.p1.set_coords(self.playground.coords(self.p1.id))

            for e in self.p1.edges:
                e.draw_VR('green')
        elif self.p2.is_active:
            self.p2.set_coords(self.playground.coords(self.p2.id))

            for e in self.p2.edges:
                e.draw_VR()

        if self.p1.tag == 'object' + str(self.obj_id):
            for v in self.p1.vertices:
                v.move(event.x - self.ex, event.y - self.ey)
        elif self.p2.tag == 'object' + str(self.obj_id):
            for v in self.p2.vertices:
                v.move(event.x - self.ex, event.y - self.ey)

        self.ex, self.ey = event.x, event.y
        self.start_v_clip()

    def start_v_clip(self):
        """
        Ready to start the v_clip algorithm
        """
        for feature in self.polygons_array[0].features:
            data = self.v_clip(self.polygons_array[0], self.polygons_array[1],
                               feature, self.polygons_array[1].features[0])

            if self.id_segment is not None:
                self.playground.delete(self.id_segment)
            if self.highlighted_feature1 is not None:
                self.playground.delete(self.highlighted_feature1)
            if self.highlighted_feature2 is not None:
                self.playground.delete(self.highlighted_feature2)

            if type(self.features_1) == Vertex and type(self.features_2) == Vertex:
                self.id_segment = self.playground.create_line(self.features_1.x, self.features_1.y,
                                                              self.features_2.x, self.features_2.y, fill='orange',
                                                              width=3)
                self.highlighted_feature1 = self.playground.create_oval(self.features_1.x - 5, self.features_1.y - 5,
                                                                        self.features_1.x + 5, self.features_1.y + 5,
                                                                        fill='yellow')

                self.highlighted_feature2 = self.playground.create_oval(self.features_2.x - 5, self.features_2.y - 5,
                                                                        self.features_2.x + 5, self.features_2.y + 5,
                                                                        fill='yellow')

            if type(self.features_1) == Vertex and type(self.features_2) == Edge:
                proj = projection(self.features_2.v2, self.features_2.v1, self.features_1)
                self.id_segment = self.playground.create_line(self.features_1.x, self.features_1.y,
                                                              proj[0], proj[1], fill='orange',
                                                              width=3)

                self.highlighted_feature1 = self.playground.create_oval(self.features_1.x - 5, self.features_1.y - 5,
                                                                        self.features_1.x + 5, self.features_1.y + 5,
                                                                        fill='yellow')

                self.highlighted_feature2 = self.playground.create_line(self.features_2.v1.x, self.features_2.v1.y,
                                                                        self.features_2.v2.x, self.features_2.v2.y,
                                                                        fill='yellow', width=3)

            if type(self.features_1) == Edge and type(self.features_2) == Edge:
                mp1 = mid_point(self.features_1.v1, self.features_1.v2)
                mp2 = mid_point(self.features_2.v1, self.features_2.v2)
                self.id_segment = self.playground.create_line(mp1, mp2, fill='orange', width=3)

                self.highlighted_feature1 = self.playground.create_line(self.features_1.v1.x, self.features_1.v1.y,
                                                                        self.features_1.v2.x, self.features_1.v2.y,
                                                                        fill='yellow', width=3)

                self.highlighted_feature2 = self.playground.create_line(self.features_2.v1.x, self.features_2.v1.y,
                                                                        self.features_2.v2.x, self.features_2.v2.y,
                                                                        fill='yellow', width=3)

    def v_clip(self, A, B, X, Y):
        """
        :param A: first polygon
        :param B: second polygon
        :param X: respective initial features X
        :param Y: respective initial features Y
        """
        self.polygon_1 = A
        self.polygon_2 = B
        self.features_1 = X
        self.features_2 = Y
        Sn = {}
        while True:
            pair = FeaturePair(self.features_1, self.features_2)
            if pair.type() == "VV":
                Sn = {FeaturePair(self.features_2, E) for E in self.polygon_2.edges if E.v1 == self.features_2 or E.v2 == self.features_2}
                if self.clip_vertex(self.features_1, self.features_2, Sn):
                    continue
                Sn = {FeaturePair(self.features_1, E) for E in self.polygon_1.edges if E.v1 == self.features_1 or E.v2 == self.features_1}
                if self.clip_vertex(self.features_2, self.features_1, Sn):
                    continue
                return [self.features_1.x - self.features_2.x, self.features_1.y - self.features_2.y]

            elif pair.type() == "VE":
                Sn = {FeaturePair(self.features_2.v1, self.features_2), FeaturePair(self.features_2.v2, self.features_2)}
                if self.clip_vertex(self.features_1, self.features_2, Sn):
                    continue
                Sn = {FeaturePair(self.features_1, E) for E in self.polygon_1.edges if E.v1 == self.features_1 or E.v2 == self.features_1}
                if self.clip_edge(self.features_2, self.features_1, Sn):
                    continue
                u = self.features_2.get_directional_vector()
                w = self.features_1 - self.features_2.v1
                q = (u[0] * w[0] + u[1] * w[1]) / (
                            u[0] * u[0] + u[1] * u[1])  # if (u[0] * u[0] + u[1] * u[1]) != 0 else 0
                c = [u[0] * q, u[1] * q]
                v = [self.features_2.v1.x + c[0], self.features_2.v1.y + c[1]]
                return [self.features_1.x - v[0], self.features_1.y - v[1]]

            elif pair.type() == "EE":
                Sn = {FeaturePair(self.features_2.v1, self.features_2), FeaturePair(self.features_2.v2, self.features_2)}
                if self.clip_edge(self.features_1, self.features_2, Sn):
                    continue

                Sn = {FeaturePair(self.features_1.v1, self.features_1), FeaturePair(self.features_1.v2, self.features_1)}
                if self.clip_edge(self.features_2, self.features_1, Sn):
                    continue

                ux, uy = self.features_1.get_directional_vector(), self.features_2.get_directional_vector()
                nx = [(uy[1] * (ux[0] * uy[1] - uy[0] * ux[1])), (ux[0] * uy[1] - uy[0] * ux[1]) * uy[0]]
                ny = [(ux[1] * (uy[0] * ux[1] - ux[0] * uy[1])), (uy[0] * ux[1] - ux[0] * uy[1]) * ux[0]]

                w = self.features_2.v1 - self.features_1.v1
                q = (ny[0] * w[0] + ny[1] * w[1]) / (ny[0] * ux[0] + ny[1] * ux[1]) if (ny[0] * ux[0] + ny[1] * ux[
                    1]) != 0 else 0
                c = [ux[0] * q, ux[1] * q]
                v = [self.features_1.v1.x + c[0], self.features_1.v1.y + c[0]]

                w1 = self.features_1.v1 - self.features_2.v1
                q1 = (nx[0] * w1[0] + nx[1] * w1[1]) / (nx[0] * uy[0] + nx[1] * uy[1]) if (nx[0] * uy[0] + nx[1] * uy[
                    1]) != 0 else 0
                c1 = [uy[0] * q1, uy[1] * q1]
                v1 = [self.features_2.v1.x + c1[0], self.features_2.v1.y + c1[0]]

                return [v[0] - v1[0], v[1] - v1[1]]

            elif pair.type() == "EV":
                self.features_1, self.features_2 = self.features_2, self.features_1
                tmp = self.polygon_1
                self.polygon_1 = self.polygon_2
                self.polygon_2 = tmp

            if self.features_2 is None:
                return None

    def clip_vertex(self, V, N, Sn):
        """
        :param V: A vertex
        :param N: A feature to be updated
        :param Sn: A set of clipping feature pairs
        :return: Test if the feature N was updated (true/false)
        """
        Sn = self.clear_all(Sn)
        for pair in Sn:
            test = self.sign_distance(self.ds(V, pair.f1, pair.f2))
            if test > 0:
                pair.f1.mark()
            else:
                pair.f2.mark()
        return self.update_clear(N, Sn)

    def ds(self, v, p1, p2):
        if p1 == p2.v1:
            p = p2.calculate_plain(p2.v1)
        else:
            p = p2.calculate_plain(p2.v2)
        return p[0][0] * v.x + p[0][1] * v.y + p[1]

    def sign_distance(self, v):
        if v <= 0:
            return -1
        return 1

    def clear_all(self, Sn):
        """
        :param Sn: Array of features
        :return: Array of updated features
        """
        for pair in Sn:
            pair.f1.clear()
            pair.f2.clear()
        return Sn

    def update_clear(self, N, Sn):
        """
        :param N: A feature to be updated
        :param Sn: A set of clipping feature pairs
        :return: Test if the feature N was updated (true/false)
        """
        M = N  # store old feature
        for pair in Sn:
            if not pair.f1.marked:
                self.update_feature(N, pair.f1)
                N = pair.f1
                break
            if not pair.f2.marked:
                self.update_feature(N, pair.f2)
                N = pair.f2
                break
        if type(M) != type(N):
            return True
        return N != M  # true if feature changed

    def update_feature(self, f, value):
        """
        :param f: Vertex/Edge
        :param value: Vertex/Edge
        :return Function updates the feature f to new feature value
        """
        try:
            if self.features_1 == f:
                self.features_1 = value
            else:
                self.features_2 = value
        except AttributeError:
            self.features_2 = value

    def clip_edge(self, E, N, Sn):
        """
        :param E: An edge
        :param N: A feature to be updated
        :param Sn: A set of clipping feature pairs
        :return: Test if the feature N was updated (true/false)
        """
        Sn = self.clear_all(Sn)
        for pair in Sn:
            d1, d2 = self.ds(E.v1, pair.f1, pair.f2), self.ds(E.v2, pair.f1, pair.f2)
            lam = d2 / (d1 - d2) if (d1 - d2) != 0 else 0
            e = E.v1 * (1 - lam) + E.v2 * lam
            u = E * self.sign_distance(d2)
            if self.sign_distance(d1 * d2) > 0:
                test = self.sign_distance(d1)
            if self.sign_distance(d1 * d2) < 0:
                v = [pair.f1.x - e[0], pair.f1.y - e[1]]
                test = self.sign_distance(u[0] * v[0] + u[1] * v[1])
            if test > 0:
                pair.f1.mark()
            else:
                pair.f2.mark()
        return self.update_clear(N, Sn)


if __name__ == '__main__':
    mw = Playground()
    mw.mainloop()
