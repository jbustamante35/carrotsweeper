import collections
import itertools
import cv2
import networkx as nx
import numpy as np
from scipy import interpolate as itp
import skimage.morphology as morphology


######
# TREE
######


class Vertex:
    def __init__(self, point, degree=0, edges=None):
        self.point = np.asarray(point)
        self.degree = degree
        self.edges = []
        self.visited = False
        if edges is not None:
            self.edges = edges

    def __str__(self):
        return str(self.point)


class Edge:
    def __init__(self, start, end=None, pixels=None):
        self.start = start
        self.end = end
        self.pixels = []
        if pixels is not None:
            self.pixels = pixels
        self.visited = False


class PosWrapper(dict):
    def __getitem__(self, key):
        return [key.point[1], key.point[0]]  # switch x and y


def build_tree(img, start=None):

    # copy image since we set visited pixels to black
    img = img.copy()
    shape = img.shape
    nWhitePixels = np.sum(img)

    # neighbor offsets (8 nbors)
    nbPxOff = np.array(
        [[-1, -1], [-1, 0], [-1, 1], [0, -1], [0, 1], [1, -1], [1, 0], [1, 1]]
    )

    queue = collections.deque()

    # a list of all graphs extracted from the skeleton
    graphs = []

    blackedPixels = 0
    # we build our graph as long as we have not blacked all white pixels!
    while nWhitePixels != blackedPixels:

        # if start not given: determine the first white pixel
        if start is None:
            it = np.nditer(img, flags=["multi_index"])
            while not it[0]:
                it.iternext()

            start = it.multi_index

        startV = Vertex(start)
        queue.append(startV)
        # print("Start vertex: ", startV)

        # set start pixel to False (visited)
        img[startV.point[0], startV.point[1]] = False
        blackedPixels += 1

        # create a new graph
        G = nx.Graph()
        G.add_node(startV)

        # build graph in a breath-first manner by adding
        # new nodes to the right and popping handled nodes to the left in queue
        while len(queue):
            currV = queue[0]
            # get current vertex
            # print("Current vertex: ", currV)

            # check all neigboor pixels
            for nbOff in nbPxOff:

                # pixel index
                pxIdx = currV.point + nbOff

                if (pxIdx[0] < 0 or pxIdx[0] >= shape[0]) or (
                    pxIdx[1] < 0 or pxIdx[1] >= shape[1]
                ):
                    continue
                    # current neigbor pixel out of image

                if img[pxIdx[0], pxIdx[1]]:
                    # print( "nb: ", pxIdx, " white ")
                    # pixel is white
                    newV = Vertex([pxIdx[0], pxIdx[1]])

                    # add edge from currV <-> newV
                    G.add_edge(currV, newV, object=Edge(currV, newV))
                    # G.add_edge(newV,currV)

                    # add node newV
                    G.add_node(newV)

                    # push vertex to queue
                    queue.append(newV)

                    # set neighbor pixel to black
                    img[pxIdx[0], pxIdx[1]] = False
                    blackedPixels += 1

            # pop currV
            queue.popleft()
        # end while

        # empty queue
        # current graph is finished ->store it
        graphs.append(G)

        # reset start
        start = None

    # end while

    return graphs, img


def get_end_nodes(g):
    return [n for n in nx.nodes(g) if nx.degree(g, n) == 1]


def merge_edges(graph):

    # copy the graph
    g = graph.copy()

    # v0 -----edge 0--- v1 ----edge 1---- v2
    #        pxL0=[]       pxL1=[]           the pixel lists
    #
    # becomes:
    #
    # v0 -----edge 0--- v1 ----edge 1---- v2
    # |_________________________________|
    #               new edge
    #    pxL = pxL0 + [v.point]  + pxL1      the resulting pixel list on the edge
    #
    # an delete the middle one
    # result:
    #
    # v0 --------- new edge ------------ v2
    #
    # where new edge contains all pixels in between!

    # start not at degree 2 nodes
    startNodes = [startN for startN in g.nodes() if nx.degree(g, startN) != 2]

    for v0 in startNodes:

        # start a line traversal from each neighbor
        startNNbs = list(nx.neighbors(g, v0))

        if not len(startNNbs):
            continue

        counter = 0
        v1 = startNNbs[counter]  # next nb of v0
        while True:

            if nx.degree(g, v1) == 2:
                # we have a node which has 2 edges = this is a line segement
                # make new edge from the two neighbors
                nbs = list(nx.neighbors(g, v1))

                # if the first neihbor is not n, make it so!
                if nbs[0] != v0:
                    nbs.reverse()

                pxL0 = g[v0][v1]["object"].pixels  # the pixel list of the edge 0
                pxL1 = g[v1][nbs[1]]["object"].pixels  # the pixel list of the edge 1

                # fuse the pixel list from right and left and add our pixel n.point
                g.add_edge(
                    v0, nbs[1], object=Edge(v0, nbs[1], pixels=pxL0 + [v1.point] + pxL1)
                )

                # delete the node n
                g.remove_node(v1)

                # set v1 to new left node
                v1 = nbs[1]

            else:
                counter += 1
                if counter == len(startNNbs):
                    break
                v1 = startNNbs[counter]  # next nb of v0

    # weight the edges according to their number of pixels
    for u, v, o in g.edges(data="object"):
        g[u][v]["weight"] = len(o.pixels)

    return g


def get_longest_path(graph, endNodes):
    """
    graph is a fully reachable graph = every node can be reached from every node
    """

    if len(endNodes) < 2:
        raise ValueError("endNodes need to contain at least 2 nodes!")

    # get all shortest paths from each endpoint to another endpoint
    allEndPointsComb = itertools.combinations(endNodes, 2)

    maxLength = 0
    maxPath = None

    for ePoints in allEndPointsComb:

        # get shortest path for these end points pairs
        sL = nx.dijkstra_path_length(graph, source=ePoints[0], target=ePoints[1])

        # dijkstra can throw if now path, but we are sure we have a path

        # store maximum
        if sL > maxLength:
            maxPath = ePoints
            maxLength = sL

    if maxPath is None:
        raise ValueError("No path found!")

    return nx.dijkstra_path(graph, source=maxPath[0], target=maxPath[1]), maxLength


def assemble_graph(binary_mask):
    """
    goes about the whole graph business
    """

    WHITE = 255

    image_bool = binary_mask == WHITE

    #################
    # SKELETONIZATION
    #################

    d = morphology.disk(2)
    img = morphology.binary_closing(image_bool, selem=d)
    skeleton = morphology.medial_axis(img)

    # Find a start pixel (not necessary)
    y = int(skeleton.shape[1] / 2)
    x = np.where(skeleton[:, y] == True)[0][0]

    # disconnect
    skeleton_copy = skeleton.copy()

    graphs, imgB = build_tree(skeleton_copy, np.array([x, y]))
    for i, g in enumerate(graphs):
        endNodes = get_end_nodes(g)
        graphs[i] = {"graph": g, "endNodes": endNodes}

    simple_graphs = []
    for g in graphs:
        newG = merge_edges(g["graph"])

        simple_graphs.append({"graph": newG, "endNodes": get_end_nodes(newG)})

    for i, g in enumerate(simple_graphs):
        path = get_longest_path(g["graph"], g["endNodes"])
        simple_graphs[i]["longestPath"] = path

    return simple_graphs


def get_shoulder_point_and_radius(image_grey):
    WHITE = 255

    ############################
    # TAKE CARE OF THAT SHOULDER
    ############################

    mask = image_grey[::-1]
    mask_trans = mask.T
    mask_reversed = mask_trans[::-1]
    shoulder_midpoint = None
    radius = None

    for i, column in enumerate(mask_reversed):
        count = collections.Counter(column)

        if count.get(WHITE, None):

            white_pixels = count.get(WHITE)
            radius = round(white_pixels / 2)

            first_white = column[::-1].argmax()

            shoulder_midpoint_y = first_white + radius
            shoulder_midpoint_x = len(mask_reversed) - i
            shoulder_midpoint = np.array([shoulder_midpoint_x, shoulder_midpoint_y])

            break

    return (shoulder_midpoint[0], shoulder_midpoint[1]), radius


def get_midline(image_grey):
    """
    determines the midline of the image
    """
    WHITE = 255

    shoulder_midpoint, shoulder_radius = get_shoulder_point_and_radius(image_grey)

    # attach a circle to the shoulder
    black_column = np.zeros((image_grey.shape[0], shoulder_radius + 10), dtype=np.uint8)
    image_grey = np.append(image_grey, black_column, axis=1)
    cv2.circle(
        image_grey,
        (shoulder_midpoint[0], shoulder_midpoint[1]),
        shoulder_radius,
        WHITE,
        -1,
    )

    simple_graphs = assemble_graph(image_grey)

    ############################
    # GET THAT MIDLINE
    ############################

    for i, g in enumerate(simple_graphs):

        graph = g["graph"]

        longestPathNodes = g["longestPath"][0]
        longestPathEdges = [
            (longestPathNodes[i], longestPathNodes[i + 1])
            for i in range(0, len(longestPathNodes) - 1)
        ]

        nodes = []
        for node in longestPathNodes:
            nodes.append(node.point)
        nodes = np.array(sorted(nodes, key=lambda x: x[1]))

        # TODO: can we do this more dynamically?
        second_last_node = nodes[-2]

        # fill in the end bit!
        # TODO: an we do this in a smarter way?
        points = np.array([second_last_node, shoulder_midpoint])
        missing_points = []
        for i in range(second_last_node[1], shoulder_midpoint[1]):
            missing_point = np.array([shoulder_midpoint[0], i])
            missing_points.append(missing_point)
        missing_points = np.array(missing_points)
        points = np.concatenate((points, missing_points))

        for e in longestPathEdges:
            px = np.array(graph[e[0]][e[1]]["object"].pixels)
            px_before_last = np.array([p for p in px if p[1] <= second_last_node[1]])
            if len(px_before_last):
                points = np.concatenate((points, px_before_last))
    points = np.array(sorted(points, key=lambda x: x[1]))

    # reverse fucking y and x to be x and y!
    points = np.array([np.array([p[1], p[0]]) for p in points])
    # return points

    # take out duplicates!
    points_set = []
    points_set_xes = []
    for i, point in enumerate(points[:, 0]):
        if point not in points_set_xes:
            points_set.append(points[i])
            points_set_xes.append(point)
    points_set = np.array(points_set)

    # return points_set

    # ########
    # # SPLINE
    # ########
    # inspiration: https://stackoverflow.com/questions/17348214/using-scipy-interpolate-splrep-function

    # x_min = points_set[:, 0].min()
    # x_max = points_set[:, 0].max()
    # new_length = x_max

    x = points_set[:, 0]
    y = points_set[:, 1]

    # new_x = np.linspace(x_min, x_max, new_length)
    new_x = np.linspace(0, 1, 1000)

    # new_y = itp.interp1d(x, y, kind="cubic")(new_x)
    # return list(zip(new_x, new_y))

    mytck, myu = itp.splprep([x, y])
    xnew, ynew = itp.splev(new_x, mytck)

    return list(zip(xnew, ynew))


def get_graph(binary_mask):
    WHITE = 255
    shoulder_midpoint, shoulder_radius = get_shoulder_point_and_radius(binary_mask)

    # attach a circle to the shoulder
    black_column = np.zeros(
        (binary_mask.shape[0], shoulder_radius + 10), dtype=np.uint8
    )
    binary_mask = np.append(binary_mask, black_column, axis=1)
    cv2.circle(
        binary_mask,
        (shoulder_midpoint[0], shoulder_midpoint[1]),
        shoulder_radius,
        WHITE,
        -1,
    )
    # return binary_mask

    simple_graphs = assemble_graph(binary_mask)

    for g in simple_graphs:

        graph = g["graph"]

        longestPathNodes = g["longestPath"][0]
        longestPathEdges = [
            (longestPathNodes[i], longestPathNodes[i + 1])
            for i in range(0, len(longestPathNodes) - 1)
        ]
        nodes = []
        for node in longestPathNodes:
            nodes.append(node.point)
        nodes = np.array(sorted(nodes, key=lambda x: x[1]))
        nodes = np.array([np.array([p[1], p[0]]) for p in nodes])

        line = None
        for e in longestPathEdges:
            px = np.array(graph[e[0]][e[1]]["object"].pixels)
            if line is None:
                line = px
            if len(px):
                line = np.concatenate((line, px))

        line = np.array(sorted(line, key=lambda x: x[1]))
        line = np.array([np.array([p[1], p[0]]) for p in line])

        return nodes, line


if __name__ == "__main__":
    #src = "/Users/creimers/Downloads/Phenotyping/carrots_new/Aal/binary_mask/{Row_9599}{Root_13}{UID_189-2018}{Genotype_B2566A}{Scale_463}{Location_California}.png"
    src = "/home/jbustamante/Dropbox/EdgarSpalding/projects/carrotsweeper/midlines/sampleimages/rawsample.png"
    mask = cv2.imread(src, cv2.IMREAD_GRAYSCALE)
    dingdong = get_graph(mask)

    #cv2.imwrite("/Users/creimers/Downloads/boingboingboing.png", dingdong)
    cv2.imwrite("/home/jbustamante/Dropbox/EdgarSpalding/projects/carrotsweeper/midlines/sampleimages/testout.png", dingdong)


