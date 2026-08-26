import numpy as np

from MJOLNIR.Geometry.GeometryConcept import  GeometryConcept, GeometryObject
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pytest

@pytest.mark.unit
def test_Concept_init():
    Concept = GeometryConcept(position=(0.0,1.0,0.0))
    assert(np.all(Concept.position==np.array([0.0,1.0,0.0])))

@pytest.mark.unit
def test_Concept_plot():
    Concept = GeometryConcept()
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    with pytest.raises(NotImplementedError):
        Concept.plot(ax)

    print(str(Concept))



# Test of GeometryObject
@pytest.mark.unit
def test_Object_init():
    GenericObject = GeometryObject(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericObject.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericObject.direction==(1.0,0.0,0.0)))


@pytest.mark.unit
def test_position():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.position = (0.0,0.0,0.0)
    print(str(GenericObject))
    assert(np.all(GenericObject.position==(0.0,0.0,0.0)))


@pytest.mark.unit
def test_Object_position_exception():
    GenericConcept = GeometryConcept(position=(0,1.0,0.0))
    with pytest.raises(AttributeError):
        GenericConcept.position=((0,0),(0,0))

    with pytest.raises(AttributeError):
        GenericConcept.position=(0,0,0,0)

@pytest.mark.unit
def test_Object_direction():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.direction = (0.0,0.0,0.5)
    assert(np.all(GenericObject.direction==(0.0,0.0,1.0)))

@pytest.mark.unit
def test_Object_direction_exception():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    with pytest.raises(AttributeError):
        GenericObject.direction=((0,0),(0,0))

    with pytest.raises(AttributeError):
        GenericObject.direction=(0,0,0)


    with pytest.raises(AttributeError):
        GenericObject.direction=(0,0,0,1)
