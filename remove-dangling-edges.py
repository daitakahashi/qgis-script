# -*- coding: utf-8 -*-

"""
***************************************************************************
* NAME: Vector network: Remove dangling edges                             *
* AUTHOR(S): Daisuke Takahashi <dtakahshi42@gmail.com>                    *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License version 3.       *
*                                                                         *
***************************************************************************
"""

from PyQt5.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterEnum,
                       QgsCoordinateReferenceSystem,
                       QgsCoordinateTransformContext,
                       QgsDistanceArea,
                       QgsGeometry)
import processing

import math
import networkx as nx
import itertools as iter

# Helper classes
class key_indexer:
    def __init__(self):
        self.counter = 0
        self.dictionary = {}
        
    def index(self, key):
        if key not in self.dictionary:
            self.counter += 1
            self.dictionary[key] = self.counter
        return self.dictionary[key]

class distance_measure_by_geom:
    def __init__(self):
        pass
    
    def measure(self, geom):
        return geom.length()
    
    def measurement_type(self):
        return 'geometry'

class distance_measure_by_ellipsoid:
    def __init__(self, source_crs):
        # Setup length measure
        self.elps_crs = QgsCoordinateReferenceSystem('EPSG:4326')
        self.trans_context = QgsCoordinateTransformContext()
        self.trans_context.calculateDatumTransforms(source_crs, self.elps_crs)
        
        self.distance_area = QgsDistanceArea()
        self.distance_area.setEllipsoid('WGS84')
        self.distance_area.setSourceCrs(source_crs, self.trans_context)
    
    def measure(self, geom):
        return self.distance_area.measureLength(geom)

    def measurement_type(self):
        return 'ellipsoid'
    
class edge_length_collecter:
    def __init__(self, distance_measure, edge_ordering):
        self.edge_list = []
        self.indexer = key_indexer()
        self.distance_measure = distance_measure
        self.edge_ordering = edge_ordering
        
    def add(self, feature):
        geom = feature.geometry()
        # only the first part will be considered
        pl = next(geom.constParts())
        geom_pl = QgsGeometry.fromPolyline(pl)
        
        index_start = self.indexer.index(pl.startPoint().asWkt())
        index_end   = self.indexer.index(pl.endPoint().asWkt())
        length      = self.distance_measure.measure(geom_pl)
        
        edge = self.edge_ordering(index_start, index_end)
        self.edge_list.append((edge[0], edge[1],
                               {'feature': feature,
                                'length': length}))
    def get(self):
        return self.edge_list
    
    def get_index_dictionary(self):
        return self.indexer.dictionary

class search_bidirectional:
    def is_dangling(self, G, v):
        return G.degree(v) == 1
    def danglings(self, G):
        return {v for (v, d) in G.degree if d == 1}
    def is_path_internode(self, G, v):
        return G.degree(v) == 2
    def inward_degree(self, G, v):
        return G.degree(v)
    def is_bidirectional(self):
        return True

class search_directional:
    def is_dangling(self, G, v):
        return G.degree(v) == 1 and G.out_degree(v) == 1
    def danglings(self, G):
        return {v for (v, d) in G.degree if d == 1 and G.out_degree(v) == 1}
    def is_path_internode(self, G, v):
        return G.in_degree(v) == 1 and G.out_degree(v) == 1
    def inward_degree(self, G, v):
        return G.out_degree(v)
    def is_bidirectional(self):
        return False

def select_all(paths, G, search_strategy):
    endpoints = {k for k in paths}
    segments  = set()
    for (ignore, ps) in paths.items():
        for x in ps:
            segments.update(x['path'])
    return (segments, endpoints)

def select_conservative(paths, G, search_strategy):
    endpoints = set()
    segments  = set()
    for (endpoint, ps) in paths.items():
        ps.sort(key=lambda x: x['L'])
        longest = ps.pop()
        if search_strategy.is_bidirectional() or search_strategy.inward_degree(G, endpoint) == 1:
            endpoints.add(longest['start'])
        for x in ps:
            segments.update(x['path'])
    return (segments, endpoints)

def edge_removal(G, threshold_length, loop_range,
                 path_selecter, search_strategy,
                 feedback):
    work = search_strategy.danglings(G)
    progress_proportion = 0
    
    for ix in loop_range:
        if feedback.isCanceled():
            raise QgsProcessingException('Cancelled.')
        
        work_count = len(work)
        work_next = set()
        dangling_paths = {}
        feedback.pushInfo('Round %d start: %d candidates' % (ix, work_count))
        
        for v in work:
            if search_strategy.is_dangling(G, v):
                adj = next(x for x in G[v])
                # This is a dangling node, so there is only one (0-th) edge
                edge_length = G[v][adj][0]['length']
            
                if edge_length < threshold_length:
                    # find total length of this branch
                    w1 = v
                    w2 = adj
                    L = 0
                
                    # Stop accumulation when it reaches to a branching point (degree > 2)
                    # or the other endpoint (degree == 1)
                    to_be_removed = set()
                    while search_strategy.is_path_internode(G, w2): # G.degree(w2) == 2:
                        to_be_removed.add((w1, w2))
                        # The node w1 is on a dangling path, and w2 is the next node
                        # (still on the dangling path).
                        # Therefore, there is only a single edge between these two nodes.
                        L += G[w1][w2][0]['length']
                        w2_old = w2
                        w2 = next(x for x in G[w2] if x != w1) ## should always exist
                        w1 = w2_old
                    # Above loop does not add the last edge. Add it.
                    to_be_removed.add((w1, w2))
                    # The node w2 is an end of the dangling path. Therefore,
                    # still there is only one edge between w1 and w2.
                    L += G[w1][w2][0]['length']
                
                    # If L < threshold, remove all edges from this leaf to the next
                    # branching point
                    if L < threshold_length:
                        if w2 in dangling_paths:
                            dangling_paths[w2].append({'L':L, 'start':v, 'path':to_be_removed})
                        else:
                            dangling_paths[w2] = [{'L':L, 'start':v, 'path':to_be_removed}]
            
        if len(dangling_paths) == 0:
            feedback.pushInfo('No dangling path found')
            ## no more edges to be removed
            break
        else:
            (r, work) = path_selecter(dangling_paths, G, search_strategy)
            if len(r) == 0:
                feedback.pushInfo('No more removable path found')
                break
            
            G.remove_edges_from(r)
            
            progress_proportion += max(0, (1 - progress_proportion)/(math.ceil(math.log2(len(r))) + 1))
            feedback.setProgress(int(progress_proportion*100) if progress_proportion < 1 else 99)
            feedback.pushInfo('Round %d: %d edges (out of %d candidates) were removed' %
                              (ix, len(r), work_count))


class RemoveDanglingEdges(QgsProcessingAlgorithm):
    """
    This algorithm takes a vector layer and removes polylines that are
      (1) dangling edges (or leaf edges),
      (2) and shorter than a given threshold length.
    In other words, it works very similar to v.clean (rmdangle).
    
    Two methods are implemented.
    1. Aggressive removal (more aggresive than v.clean)
    This method removes all dangling edges in parallel (note, v.clean removes
    them sequentially), and iterate the removal step until there are no more
    dangling edges shorter than the threshold.
    
    After this method, sub-trees consisting of short dangling edges will be
    removed completely. Let the threshold be 6 and consider below subgraph
    (condsider '-' being a unit length dangling edges and '=' be
    non-dangling edges)
    
                      O --- O
                     /
        ~ === O === O - O
               \
                O
    
    After a single removal step, we will get below topology. These nodes and edges
    may be removed in the following iterations.
    
        ~ --- O --- O
    
    For each iteration, the number of candidate nodes will reduce exponentially.
    Therefore, this method does not need many iterations and is quite efficient
    (but may remove too much detail).
    
    2. Conservative removal
    This method takes two pass approach to keep a long subpath inside of a sub-tree
    of short dangling edges. Let the threshold be 6, and consider the same subgraph
    as above.
    
                     (*)---(*)
                     /
        ~ ===(*)===(*)- O
               \
                O
    
    In this subgraph, a path connecting nodes (*) is longer than
    the threshold 6, thus, and this method will keep it.
    
    At the first path, this method removes dangling paths only from nodes having
    exactly one non-dangling edge (edge 1 and 2 in the below figure).
    
                      O --- O
                   1 /
        ~ === O === O - O
             3 \      2
                O
    
    Then, group these dangling paths by their endpoints, and for each group,
    remove all dangling paths except the longest one.
    
                      O --- O
                   1 /
        ~ === O --- O 
             3 \
                O
    
    Now the path 1 and 3 join at a node that has exactly one non-dangling edge,
    thus, these paths will be considered in the next iteration.
    The method iterate above removal step until no more edges are to be removed.
    
    After the first pass, we will get a graph like below.
   
    ~ === O === O === O === ~
           \     \     
            \     O
             \
              O --- O --- O
    
    These 'wiskers' are then removed by applying "aggresive' method once (pass 2).
    
    In the worst case scenario, the pass 1 reduces the number of candidates only
    linearly, and it may take a lot of iterations to complete the transformation.
    
    A result of this conservative method will be the same as iterative applications
    of v.clean with gradually increasing thresholds (with very small steps).
    
    Note: Both methods do not consider edge direction.
    Loops will never be removed.
    """
    
    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT = 'INPUT'
    THRESHOLD_LENGTH = 'THRESHOLD_LENGTH'
    USE_ELLIPSOID = 'USE_ELLIPSOID'
    KEEP_LONG = 'KEEP_LONG'
    DIRECTION = 'DIRECTION'
    MAX_ITER = 'MAX_ITER'
    
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return RemoveDanglingEdges()

    def name(self):
        return 'removedanglingedges'

    def displayName(self):
        return self.tr('Remove dangling edges')

    def group(self):
        return self.tr('Vector')

    def groupId(self):
        return 'vector'

    def shortHelpString(self):
        return self.tr("Iteratively remove dangling edges")

    def initAlgorithm(self, config=None):
        # We add the input vector features source. It have to be a polyline layer.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorLine]
            )
        )
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.THRESHOLD_LENGTH,
                self.tr('Threhsold length of a dangling edge to be removed'),
                QgsProcessingParameterNumber.Double,
                0,
                False,
                0
            )
        )
        
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.USE_ELLIPSOID,
                self.tr('Use ellipsoid for length calculation (otherwise use a map unit)'),
                False
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.KEEP_LONG,
                self.tr('Merge dangling edges first to keep long paths inside'),
                True
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.DIRECTION,
                self.tr('Search direction'),
                ['bidirectional', 'from start', 'from end'],
                defaultValue=0
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_ITER,
                self.tr('Maximum number of removal iterations (0 to be unlimited)'),
                QgsProcessingParameterNumber.Integer,
                0,
                False,
                0
            )
        )

        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output layer'),
                QgsProcessing.TypeVectorLine
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        # Retreave source and sink layers
        source = self.parameterAsSource(
            parameters,
            self.INPUT,
            context
        )
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            source.fields(),
            source.wkbType(),
            source.sourceCrs()
        )
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        # Retreave other parameters
        threshold_length = self.parameterAsDouble(
            parameters,
            self.THRESHOLD_LENGTH,
            context
        )

        use_ellipsoid = self.parameterAsBool(
            parameters,
            self.USE_ELLIPSOID,
            context
        )
        
        keep_long = self.parameterAsBool(
            parameters,
            self.KEEP_LONG,
            context
        )

        direction = self.parameterAsEnum(
            parameters,
            self.DIRECTION,
            context
        )
        
        maximum_iteration = self.parameterAsInt(
            parameters,
            self.MAX_ITER,
            context
        )
        
        fc = source.featureCount()
        total = 100.0 / source.featureCount() if fc > 0 else 0
        feedback.pushInfo('%d features found' % fc)
        
        feedback.setProgressText('Building a graph...')
        feedback.setProgressText('Collecting edges...')
        
        if use_ellipsoid:
            distance_measure = distance_measure_by_ellipsoid(source.sourceCrs())
        else:
            distance_measure = distance_measure_by_geom()

            
        # G sholud be a multigraph, which is a graph allowing multiple edges
        # between two nodes, because several lines may connect the same endpoints
        # through different ways.
        #
        # Choose appropriate graph type and ordering according to the specified
        # option.
        if direction == 0:
            G = nx.MultiGraph()
            search_strategy = search_bidirectional()
        else:
            G = nx.MultiDiGraph()
            search_strategy = search_directional()
        
        # When finding danglings from endpoints, build a graph with reverseing edges
        # directions.
        if direction == 2:
            edge_collecter = edge_length_collecter(distance_measure,
                                                   lambda x, y: (y, x))
        else:
            edge_collecter = edge_length_collecter(distance_measure,
                                                   lambda x, y: (x, y))
        
        features = source.getFeatures()
        for (ix, f) in enumerate(features):
            if feedback.isCanceled():
                raise QgsProcessingException('Cancelled.')
                
            edge_collecter.add(f)
            feedback.setProgress(int(ix * total))
            
        feedback.setProgressText('Add edges to a graph...')
            
        G.add_edges_from(edge_collecter.get())
        feedback.pushInfo('The graph has %d edges and %d nodes' %
                          (G.number_of_edges(), G.number_of_nodes()))
        
        feedback.setProgress(0)

        feedback.pushInfo('Dangling edges shorter than %d (using %s) will be removed' %
                          (threshold_length, distance_measure.measurement_type))
        
        feedback.setProgressText('Removing dangling edges')
        
        if maximum_iteration == 0:
            loop_range = iter.count(1)
        else:
            loop_range = range(1, maximum_iteration + 1)
        
        # Main process
        if keep_long:
            feedback.pushInfo('Conservative removal: Keep subpaths longer than the threshold')
            feedback.setProgressText('Pass 1: Merge trees of short dangling edges')
            edge_removal(G, threshold_length, loop_range,
                         select_conservative, search_strategy,
                         feedback)
            feedback.setProgressText('Pass 2: Remove remaining dangling edges (only once)')
            edge_removal(G, threshold_length, [1],
                         select_all, search_strategy,
                         feedback)
        else:
            feedback.pushInfo('Agressive removal: Remove all short dangling edges for every round')
            edge_removal(G, threshold_length, loop_range,
                         select_all, search_strategy,
                         feedback)
        
        feedback.setProgressText('Copying %d remaining features...' % G.number_of_edges())
        features = [G.edges[e]['feature'] for e in G.edges]
        sink.addFeatures(features, QgsFeatureSink.FastInsert)
        
        feedback.setProgressText('Done.')
        
        return {self.OUTPUT: dest_id}

