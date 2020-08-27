import numpy as np

from MJOLNIR.Geometry.GeometryConcept import  GeometryConcept, GeometryObject
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def test_Concept_init():
    Concept = GeometryConcept(position=(0.0,1.0,0.0))
    assert(np.all(Concept.position==np.array([0.0,1.0,0.0])))


def test_Concept_plot():
    Concept = GeometryConcept()
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    try:
        Concept.plot(ax)
        assert False
    except NotImplementedError:
        assert True
    print(str(Concept))



# Test of GeometryObject
def test_Object_init():
    GenericObject = GeometryObject(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericObject.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericObject.direction==(1.0,0.0,0.0)))

def test_position():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.position = (0.0,0.0,0.0)
    print(str(GenericObject))
    assert(np.all(GenericObject.position==(0.0,0.0,0.0)))

def test_Object_position_exception():
    GenericConcept = GeometryConcept(position=(0,1.0,0.0))
    try:
        GenericConcept.position=((0,0),(0,0))
        assert False
    except AttributeError:
        assert True

    try:
        GenericConcept.position=(0,0,0,0)
        assert False
    except AttributeError:
        assert True

def test_Object_direction():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.direction = (0.0,0.0,0.5)
    assert(np.all(GenericObject.direction==(0.0,0.0,1.0)))

def test_Object_direction_exception():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    try:
        GenericObject.direction=((0,0),(0,0))
        assert False
    except AttributeError:
        assert True

    try:
        GenericObject.direction=(0,0,0)
        assert False
    except AttributeError:
        assert True

    try:
        GenericObject.direction=(0,0,0,1)
        assert False
    except AttributeError:
        assert True