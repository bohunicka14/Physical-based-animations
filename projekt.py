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

class VoronoiRegion:

    def __init__(self, canvas, vertex, endpoint):
        self.canvas = canvas
        self.vertex = vertex
        self.endpoint = endpoint

    def draw(self):
        self.canvas.create_line(self.vertex + self.endpoint, fill = 'red')


class Feature():

    def mark(self):
        self.marked = True

    def clear(self):
        self.marked = False

class Vertex(Feature):

    def __init__(self, x, y, marked = False):
        self.x = x
        self.y = y
        self.marked = marked
        self.voronoi_region = []

    def draw(self, canvas):
        self.canvas = canvas
        self.id = self.canvas.create_oval(self.x-self.r, self.y-self.r, self.x+self.r, self.y+self.r,
                                fill='black')

    def update(self):
        self.canvas.itemconfig(self.id, fill='red')

    def move(self, xx, yy):
        self.canvas.move(self.id, xx, yy)

class Edge(Feature):

    def __init__(self, v1, v2, marked = False):
        self.v1 = v1
        self.v2 = v2
        self.marked = marked
        self.voronoi_region = []
        
    def update(self, canvas):
        self.canvas = canvas
        self.canvas.create_line(self.v1.x, self.v1.y, self.v2.x, self.v2.y, fill='red')

    # def update(self, x1, y1, x2, y2):
    #     self.v1.update(x1, y1)
    #     self.v2.update(x2, y2)

    def get_directional_vector(self):
        return [self.v2.x - self.v1.x, self.v2.y - self.v1.y]

    def get_normal_vector(self):
        directional = self.get_directional_vector()
        return [directional[1], -directional[0]]

    def get_inverse_normal_vector(self):
        normal = self.get_normal_vector()
        return [-normal[0], -normal[1]]

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
        self.features = []

    def draw(self, coords, color):
        """ Coords - budem posielat suradnice bodov na vykreslenie poly-lineu
        """
        # print(type(coords))
        if self.is_active:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='red', width=3, tag=self.tag)
        else:
            self.id = self.canvas.create_polygon(*coords, fill=color, outline='black', width=1, tag=self.tag)

    def draw_vertices(self):
        for vertex in self.vertices:
            vertex.draw(self.canvas)

    def set_features(self):
        if len(self.vertices) == 0 and len(self.edges) == 0:
            ## init features
            for i in range(0, len(self.coords), 2):
                # print(len(self.coords), self.coords[i], self.coords[i+1])
                self.vertices.append(Vertex(self.coords[i], self.coords[i + 1]))
            self.draw_vertices()
            for i in range(len(self.vertices) - 1):
                self.edges.append(Edge(self.vertices[i], self.vertices[i + 1]))
            self.edges.append(Edge(self.vertices[0], self.vertices[-1]))
            # print(len(self.edges))
        else:
            ## update features according to new coords
            j = 0
            for i in range(0, len(self.coords), 2):
                self.vertices[j].update(self.coords[i], self.coords[i+1])
                j += 1
        self.features = self.vertices + self.edges

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
        # print('Suradnice bodov objektu p1: ', self.playground.coords(self.p1.id))
        self.p2 = PolyObject(self.playground, 400, 400, 'object2')
        self.polygons_array.append(self.p1)
        self.polygons_array.append(self.p2)
        # print(self.polygons_array)

        self.p1.draw(coords_1, 'white')
        self.p2.draw(coords_2, 'white')
        self.p1.set_coords(self.playground.coords(self.p1.id))
        self.p2.set_coords(self.playground.coords(self.p2.id))
        print('coords   ', self.p1.coords)

        self.v_clip(self.p1, self.p2, None, None)

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
        self.p1.get_centroid()
        self.p2.get_centroid()

    def drag(self, event):
        # if self.obj_tag == ('object' + str(self.obj_id)):
        self.playground.move(self.obj_id, event.x - self.ex, event.y - self.ey)
        if self.p1.tag == 'object' + str(self.obj_id):
            for v in self.p1.vertices:
                v.move(event.x - self.ex, event.y - self.ey)
        elif self.p2.tag == 'object' + str(self.obj_id):
            for v in self.p2.vertices:
                v.move(event.x - self.ex, event.y - self.ey)

        # if self.p1.is_active:
            # self.p1.set_coords(self.playground.coords(self.p1.id))
        # elif self.p2.is_active:
            # self.p2.set_coords(self.playground.coords(self.p2.id))
        self.start_v_clip()
        self.ex, self.ey = event.x, event.y
        
    def start_v_clip(self):
        """
        Ready to start the v_clip algorithm
        """
        # for feature in self.polygons_array[0].features:
        #     # print('som tu', self.polygons_array[0], self.polygons_array[1],
        #     #                    feature, self.polygons_array[1].features)
        #     self.v_clip(self.polygons_array[0], self.polygons_array[1],
        #                        feature, self.polygons_array[1].features[0])

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
        """
        :param V: A vertex
        :param N: A feature to be updated
        :param Sn: A set of clipping feature pairs
        :return: Test if the feature N was updated (true/false)
        """
        self.clear_all(Sn)
        for pair in Sn:
            test = self.sign_distance(self.ds(V, self.voronoi_plane(pair[0], pair[1])))
            if test > 0:
                pair.f1.mark()
            else:
                pair.f2.mark()
        return self.update_clear(N, Sn)

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
        M = N                   # store old feature
        for pair in Sn:
            if not pair.f1.marked:
                N = pair.f1
            if not pair.f2.marked:
                N = pair.f2
        return N != M           # true if feature changed

    def clip_edge(self, E, N, Sn):
        """
        :param E: An edge
        :param N: A feature to be updated
        :param Sn: A set of clipping feature pairs
        :return: Test if the feature N was updated (true/false)
        """
        self.clear_all(Sn)
        for x, y in Sn:
            pass
        return self.update_clear(N, Sn)


if __name__ == '__main__':
    mw = Playground()
    mw.mainloop()
