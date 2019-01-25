import tkinter
from tkinter import messagebox
from tkinter import *
import math


def multiply(m, array):
    for i in array:
        i = int(i)

    for i in range(0, len(array), 2):
        result_x = array[i] * m[0][0] + array[i + 1] * m[0][1]
        result_y = array[i] * m[1][0] + array[i + 1] * m[1][1]
        array[i] = result_x
        array[i + 1] = result_y
    return array

def unit_vector(vector):
    length = math.sqrt(vector[0]**2 + vector[1]**2)
    return [vector[0] / length, vector[1] / length]

def ray_line_segment_intersection(o, d, a, b):
    ## returns true if ray intersects with line segment
    ## o - point of ray
    ## d - dir vector of ray
    ## a,b - line segment points
    ortho1 = [-d[1], d[0]]
    ortho2 = [d[1], -d[0]]
    a_to_o = [o[0] - a[0], o[1] - a[1]]
    a_to_b = [b[0] - a[0], b[1] - a[1]]

    denom1 = a_to_b[0] * ortho1[0] + a_to_b[1] * ortho1[1]
    if denom1 == 0:
        return False

    denom2 = a_to_b[0] * ortho2[0] + a_to_b[1] * ortho2[1]
    if denom2 == 0:
        return False

    t1 = abs(a_to_b[0] * a_to_o[1] - a_to_o[0] * a_to_b[1]) / denom1
    t2 = (a_to_o[0] * ortho1[0] + a_to_o[1] * ortho1[1]) / denom1

    t12 = abs(a_to_b[0] * a_to_o[1] - a_to_o[0] * a_to_b[1]) / denom2
    t22 = (a_to_o[0] * ortho2[0] + a_to_o[1] * ortho1[1]) / denom2

    partial_result1 = t2 >= 0 and t2 <= 1 and t1 >= 0
    partial_result2 = t22 >= 0 and t22 <= 1 and t12 >= 0


    return partial_result1 and partial_result2


# class VoronoiRegion:
#
#     def __init__(self, canvas, vertex1, endpoint1, vertex2, endpoint2):
#         self.canvas = canvas
#         self.vertex1 = vertex1
#         self.endpoint1 = endpoint1
#         self.vertex2 = vertex2
#         self.endpoint2 = endpoint2
#
#     def draw(self):
#         self.canvas.create_line(self.vertex1 + self.endpoint1, fill = 'red')
#         self.canvas.create_line(self.vertex2 + self.endpoint2, fill = 'red')


class Feature():

    def mark(self):
        self.marked = True

    def clear(self):
        self.marked = False

class Vertex(Feature):

    def __init__(self, canvas, x, y, marked = False):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.marked = marked
        self.voronoi_region = [[self.x, self.y], [self.x, self.y]]
        self.vr1 = self.vr2 = None

    def update(self, x, y):
        self.x = x
        self.y = y
        self.voronoi_region = [[self.x, self.y], [self.x, self.y]]

    def draw_VR(self, color = 'blue'):
        if self.vr2 is not None and self.vr1 is not None:
            self.canvas.delete(self.vr1)
            self.canvas.delete(self.vr2)
        if len(self.voronoi_region[0]) == 4:
            self.vr1 = self.canvas.create_line(self.voronoi_region[0], fill=color)
        if len(self.voronoi_region[1]) == 4:
            self.vr2 = self.canvas.create_line(self.voronoi_region[1], fill=color)

    def move_vr(self, x, y):
        self.canvas.move(self.vr1, x, y)
        self.canvas.move(self.vr2, x, y)


class Edge(Feature):

    def __init__(self, v1, v2, polygon, marked = False):
        self.v1 = v1
        self.v2 = v2
        self.polygon = polygon
        self.marked = marked
        self.voronoi_region = self.create_voronoi_region()
        self.vr1 = self.vr2 = None

    # def update(self, x1, y1, x2, y2):
    #     self.v1.update(x1, y1)
    #     self.v2.update(x2, y2)

    def get_directional_vector(self):
        return [self.v2.x - self.v1.x, self.v2.y - self.v1.y]

    def get_normal_vector(self):
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
        normal = self.get_normal_vector()
        return [-normal[0], -normal[1]]

    def create_voronoi_region(self):
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

    def draw_VR(self, color = 'blue'):
        if self.vr2 is not None and self.vr1 is not None:
            self.polygon.canvas.delete(self.vr1)
            self.polygon.canvas.delete(self.vr2)
        self.vr1 = self.polygon.canvas.create_line(self.voronoi_region[0], fill=color)
        self.vr2 = self.polygon.canvas.create_line(self.voronoi_region[1], fill=color)

    def move_vr(self, x, y):
        self.polygon.canvas.move(self.vr1, x, y)
        self.polygon.canvas.move(self.vr2, x, y)

    def update(self):
        self.voronoi_region = self.create_voronoi_region()


class FeaturePair():

    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2

    def type(self):
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
        tmp = self.f1
        self.f1 = self.f2
        self.f2 = tmp



class PolyObject:
    ### !!!!!!! Pridat pole alebo zoznam so zapamatanymi hranami, ktore spajaju vrcholy a pamatat si vrcholy
    def __init__(self, canvas, x, y, tag):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.tag = tag
        self.is_active = 0
        self.edges = []
        self.vertices = []
        self.voronoi_regions = []

    def draw(self, coords, color):
        """ Coords - budem posielat suradnice bodov na vykreslenie poly-lineu
        """
        # print(type(coords))
        if self.is_active:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='red', width=3, tag=self.tag)
        else:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='black', width=1, tag=self.tag)

    def point_in_polygon(self, x, y):
        n = len(self.vertices)
        v = False

        b1x, b1y = self.vertices[0].x, self.vertices[0].y
        for i in range(n+1):
            b2x, b2y = self.vertices[i%n].x, self.vertices[i%n].y
            if y > min(b1y, b2y):
                if y <= max(b1y, b2y):
                    if x <= max(b1x, b2x):
                        if b1y != b2y:
                            xin = (y-b1y)*(b2x-b1x)/(b2y-b1y) + b1x
                        if b1x == b2x or x <= xin:
                            v = not v
            b1x,b1y = b2x,b2y
        return v

    def set_active(self, state):
        """ Zmena outline farby
        """
        if state == 1:
            self.canvas.itemconfig(self.id, outline='red', width=3)
            self.is_active = 1
        else:
            self.canvas.itemconfig(self.id, outline='black', width=1)
            self.is_active = 0

    def set_features(self):
        if len(self.vertices) == 0 and len(self.edges) == 0:
            ## init features
            for i in range(0, len(self.coords), 2):
                # print(len(self.coords), self.coords[i], self.coords[i+1])
                self.vertices.append(Vertex(self.canvas, self.coords[i], self.coords[i + 1]))
            for i in range(len(self.vertices) - 1):
                self.edges.append(Edge(self.vertices[i], self.vertices[i + 1], self))
            self.edges.append(Edge(self.vertices[0], self.vertices[-1], self))
            # print(len(self.edges))
        else:
            ## update features according to new coords
            j = 0
            for i in range(0, len(self.coords), 2):
                self.vertices[j].update(self.coords[i], self.coords[i+1])
                j += 1

            for e in self.edges:
                e.update()

    def set_coords(self, coords):
        self.coords = coords
        self.set_features()

    def get_centroid(self):
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
        # matrix = [[1, 0, 0],
        # 		  [0, 1, 0],
        # 		  [xx - self.x, yy - self.y, 1]]
        # self.x = xx - self.x
        # self.y = yy - self.y
        # self.coords = multiply(matrix, self.coords)
        self.canvas.move(self.id, xx, yy)


class Playground(Tk):
    def __init__(self):
        super().__init__()
        main = Frame(self, bg='red', width=1000)
        main.pack()

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

        if self.p1.is_active:
            # print(self.playground.coords(self.p1.id))
            self.obj_id = self.p1.id
            self.p1.rotate()
        elif self.p2.is_active:
            # print(self.playground.coords(self.p2.id))
            self.obj_id = self.p2.id
            # self.p2.set_coords(rotated_points)
            self.p2.rotate()
        else:
            messagebox.showwarning("Warning", "No object is selected!!!")
            return

        self.redraw_canvas()

    def load_from_file(self, file):
        """ Naciatnie suradnic zo suboru
        """
        coords = []
        f = open(file, 'r')

        for line in f:
            coords.append([line.strip(' ')])
        f.close()
        self.create_polygons(coords)

    def create_polygons(self, polygons):
        """ Z nacitanych suradnic vytvor objekty polygonov a vykresli ich
        """
        coords_1 = polygons[0]
        coords_2 = polygons[1]
        self.p1 = PolyObject(self.playground, 100, 100, 'object1')
        self.p1.draw(coords_1, 'white')
        # print('Suradnice bodov objektu p1: ', self.playground.coords(self.p1.id))
        self.p2 = PolyObject(self.playground, 400, 400, 'object2')
        self.p2.draw(coords_2, 'navy')
        self.polygons_array.append(self.p1)
        self.polygons_array.append(self.p2)
        # print(self.polygons_array)

        self.p1.set_coords(self.playground.coords(self.p1.id))
        self.p2.set_coords(self.playground.coords(self.p2.id))

        # # for v in self.p1.vertices:
        # #     v.draw_VR()
        #
        # for e in self.p1.edges:
        #     e.draw_VR()
        #
        # # for v in self.p2.vertices:
        # #     v.draw_VR()
        #
        # for e in self.p2.edges:
        #     e.draw_VR()


        # print('coords   ', self.p1.coords)

        # self.v_clip(self.p1, self.p2, None, None)

    def redraw_canvas(self):
        self.playground.delete('all')
        self.p1.draw(self.p1.coords, 'white')
        self.p2.draw(self.p2.coords, 'navy')


    def click(self, event):
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

        # self.initial_coords = self.playground.coords(self.obj_id)
        self.ex, self.ey = event.x, event.y

    def drop(self, event):
        if self.p1.is_active:
            self.p1.set_coords(self.playground.coords(self.p1.id))
        elif self.p2.is_active:
            self.p2.set_coords(self.playground.coords(self.p2.id))
        # self.p1.get_centroid()
        # self.p2.get_centroid()

    def drag(self, event):
        # if self.obj_tag == ('object' + str(self.obj_id)):
        self.playground.move(self.obj_id, event.x - self.ex, event.y - self.ey)
        if self.p1.is_active:
            self.p1.set_coords(self.playground.coords(self.p1.id))
            # for e in self.p1.edges:
            #     e.move_vr(event.x - self.ex, event.y - self.ey)
            for v in self.p1.vertices:
                v.draw_VR()

            for e in self.p1.edges:
                e.draw_VR()
        elif self.p2.is_active:
            self.p2.set_coords(self.playground.coords(self.p2.id))
            # for e in self.p2.edges:
            #     e.move_vr(event.x - self.ex, event.y - self.ey)
            for v in self.p2.vertices:
                v.draw_VR()

            for e in self.p2.edges:
                e.draw_VR()

        self.ex, self.ey = event.x, event.y


    def v_clip(self, A, B, X, Y):
        # print(A, B, X, Y)
        pair = FeaturePair(X, Y)
        Sn = {}
        while True:
            if pair.type() == "VV":
                Sn = {FeaturePair(Y, E) for E in B.edges if E.v1 == Y or E.v2 == Y}
                if self.clip_vertex(X, Y, Sn):
                    continue
                Sn = {FeaturePair(X, E) for E in A.edges if E.v1 == X or E.v2 == X}
                if self.clip_vertex(Y, X, Sn):
                    continue
                return [X.x - Y.x, X.y - Y.y]

            elif pair.type() == "VE":
                Sn = {FeaturePair(Y.v1, Y), FeaturePair(Y.v2, Y)}
                if self.clip_vertex(X, Y, Sn):
                    continue
                Sn = {FeaturePair(X, E) for E in A.edges if E.v1 == X or E.v2 == X}
                if self.clip_edge(Y, X, Sn):
                    continue
                u = [Y.v2.x - Y.v1.x, Y.v2.y - Y.v1.y]

                return [X.x - (Y.v1.x + ((u[0] * (X.x - Y.v1.x)) / (u[0] * u[0])) * u[0]),
                        X.y - (Y.v1.y + ((u[1] * (X.y - Y.v1.y)) / (u[1] * u[1])) * u[1])]

            elif pair.type() == "EE":
                Sn = {FeaturePair(Y.v1, Y), FeaturePair(Y.v2, Y)}
                if self.clip_edge(X, Y, Sn):
                    continue

                Sn = {FeaturePair(X.v1, X), FeaturePair(X.v2, X)}
                if self.clip_edge(Y, X, Sn):
                    continue

                ux = [X.v2.x - X.v1.x, X.v2.y - X.v1.y]
                uy = [Y.v2.x - Y.v1.x, Y.v2.y - Y.v1.y]
                nx = [ux[1], -ux[0]]
                ny = [uy[1], -uy[0]]

                tmp1 = [X.v1.x + ((ny[0] * (Y.v1.x - X.v1.x)) / (ny[0] * ux[0])) * ux[0],
                        X.v1.y + ((ny[1] * (Y.v1.y - X.v1.y)) / (ny[1] * ux[1])) * ux[1]]

                tmp2 = [Y.v1.x + ((nx[0] * (X.v1.x - Y.v1.x)) / (nx[0] * uy[0])) * uy[0],
                        Y.v1.y + ((nx[1] * (X.v1.y - Y.v1.y)) / (nx[1] * uy[1])) * uy[1]]
                return [tmp1[0] - tmp2[0], tmp1[1] - tmp2[1]]

            elif pair.type() == "EV":
                pair.swap() # swap(X, Y)
                # swap(A, B)
                tmp = A
                A = B
                B = tmp

            if Y is None:
                return None


    def clip_vertex(self, V, N, Sn):
        # V - vrchol, N - cast, ktora bude updatovana, Sn - set of clipping feature pairs
        # return bool , cize true/false
        self.clear_all(Sn)
        for x, y in Sn:
            pass
        return self.update_clear(N, Sn)

    def update_clear(self, N, Sn):
        # return bool - test if the feature N was updated
        M = N                   # store old feature
        for x, y in Sn:
            pass
        return N != M           # true if feature changed

    def clip_edge(self, E, N, Sn):
        self.clear_all(Sn)
        for x, y in Sn:
            pass
        return self.update_clear(N, Sn)


if __name__ == '__main__':
    mw = Playground()
    mw.mainloop()
