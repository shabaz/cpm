#include <Python.h>
#include <numpy/arrayobject.h>

#include "cpm.h"

typedef struct {
    PyObject_HEAD
    Cpm<Lattice2d>* ptrObj;
} PyCpm2d;

typedef struct {
    PyObject_HEAD
    Cpm<Lattice3d>* ptrObj;
} PyCpm3d;



static PyModuleDef cpmmodule = {
    PyModuleDef_HEAD_INIT,
    "cpm",
    "CPM Python wrapper",
    -1,
    NULL, NULL, NULL, NULL, NULL
};


static int PyCpm2d_init(PyCpm2d *self, PyObject* args, PyObject* kwds) {
    int dimension;
    int numberOfTypes;
    int temperature;
    if (! PyArg_ParseTuple(args, "iii", &dimension, &numberOfTypes, &temperature))
        return -1;

    self->ptrObj = new Cpm<Lattice2d>(dimension, numberOfTypes, temperature);
    return 0;
}

static void PyCpm2d_dealloc(PyCpm2d* self) {
    delete self->ptrObj;
    Py_TYPE(self)->tp_free(self);
}



static PyTypeObject PyCpm2dType = { PyVarObject_HEAD_INIT(NULL, 0)
                                    "cpm.Cpm2d"   /* tp_name */
                                };


static int PyCpm3d_init(PyCpm3d *self, PyObject* args, PyObject* kwds) {
    int dimension;
    int numberOfTypes;
    int temperature;
    if (! PyArg_ParseTuple(args, "iii", &dimension, &numberOfTypes, &temperature))
        return -1;

    self->ptrObj = new Cpm<Lattice3d>(dimension, numberOfTypes, temperature);
    return 0;
}

static void PyCpm3d_dealloc(PyCpm3d* self) {
    delete self->ptrObj;
    Py_TYPE(self)->tp_free(self);
}



static PyTypeObject PyCpm3dType = { PyVarObject_HEAD_INIT(NULL, 0)
                                    "cpm.Cpm3d"   /* tp_name */
                                };

static PyObject * PyCpm2d_setConstraints(PyCpm2d* self, PyObject* args, 
        PyObject* kwargs )
{
    char* keywords [] = {
        "cell_type", 
        "other_cell_type", 
        "lambda_perimeter", 
        "target_perimeter",
        "lambda_area", 
        "target_area",
        "lambda_act",
        "max_act",
        "adhesion",
        "lambda_connectedness",
        "lambda_persistence",
        "persistence_diffusion",
        "persistence_time",
        "fixed",
        "lambda_chemotaxis",
        NULL
    };
    int cellType = -1, targetPerimeter = -1, 
        targetArea = -1, maxAct = -1, 
        otherCellType = -1, adhesion = -1, persistenceTime = -1, fixed = -1;
    double lambdaPerimeter = -1, lambdaArea = -1, lambdaAct = -1, connectedLambda = -1, 
           persistenceLambda = -1, persistenceDiffusion = -1, lambdaChemotaxis=-1;
     if (! PyArg_ParseTupleAndKeywords(args, kwargs, "i|$idididiidddiid", 
                 keywords, &cellType, &otherCellType, &lambdaPerimeter, 
                 &targetPerimeter, &lambdaArea, &targetArea, &lambdaAct, 
                 &maxAct, &adhesion, &connectedLambda, &persistenceLambda, 
                 &persistenceDiffusion, &persistenceTime, &fixed, &lambdaChemotaxis))
         return Py_False;

     if (fixed >= 0) {
         (self->ptrObj)->setFixedConstraint(cellType, fixed);
     }

     if (lambdaPerimeter >= 0 && targetPerimeter != -1) {
         (self->ptrObj)->setPerimeterConstraints(cellType, lambdaPerimeter, 
                 targetPerimeter);
     }

     if (lambdaArea >= 0 && targetArea != -1) {
         (self->ptrObj)->setAreaConstraints(cellType, lambdaArea, targetArea);
     }

     if (lambdaAct >= 0 && maxAct != -1) {
         (self->ptrObj)->setActConstraints(cellType, lambdaAct, maxAct);
     }

     if (connectedLambda >= 0) {
         (self->ptrObj)->setConnectedConstraints(cellType, connectedLambda);
     }

     if (lambdaChemotaxis >= 0) {
         (self->ptrObj)->setChemotaxisConstraints(cellType, lambdaChemotaxis);
     }

     if (persistenceLambda >= 0 && persistenceDiffusion >= 0 && persistenceTime >= 0) {
         (self->ptrObj)->setPersistenceConstraints(cellType, persistenceLambda, 
                 persistenceTime, persistenceDiffusion);
     }

     if (otherCellType != -1 && adhesion != -1) {
         (self->ptrObj)->setAdhesionBetweenTypes(cellType, otherCellType, 
                 adhesion);
     }

    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject * PyCpm3d_setConstraints(PyCpm3d* self, PyObject* args, 
        PyObject* kwargs )
{
    char* keywords [] = {
        "cell_type", 
        "other_cell_type", 
        "lambda_perimeter", 
        "target_perimeter",
        "lambda_area", 
        "target_area",
        "lambda_act",
        "max_act",
        "adhesion",
        "lambda_connectedness",
        "lambda_persistence",
        "persistence_diffusion",
        "persistence_time",
        "fixed",
        NULL
    };
    int cellType = -1, targetPerimeter = -1, 
        targetArea = -1, maxAct = -1, 
        otherCellType = -1, adhesion = -1, persistenceTime = -1, fixed = -1;
    double lambdaPerimeter = -1, lambdaArea = -1, lambdaAct = -1, connectedLambda = -1, 
           persistenceLambda = -1, persistenceDiffusion = -1;
     if (! PyArg_ParseTupleAndKeywords(args, kwargs, "i|$idididiidddii", 
                 keywords, &cellType, &otherCellType, &lambdaPerimeter, 
                 &targetPerimeter, &lambdaArea, &targetArea, &lambdaAct, 
                 &maxAct, &adhesion, &connectedLambda, &persistenceLambda, 
                 &persistenceDiffusion, &persistenceTime, &fixed))
         return Py_False;

     if (fixed >= 0) {
         (self->ptrObj)->setFixedConstraint(cellType, fixed);
     }

     if (lambdaPerimeter >= 0 && targetPerimeter != -1) {
         (self->ptrObj)->setPerimeterConstraints(cellType, lambdaPerimeter, 
                 targetPerimeter);
     }

     if (lambdaArea >= 0 && targetArea != -1) {
         (self->ptrObj)->setAreaConstraints(cellType, lambdaArea, targetArea);
     }

     if (lambdaAct >= 0 && maxAct != -1) {
         (self->ptrObj)->setActConstraints(cellType, lambdaAct, maxAct);
     }

     if (connectedLambda >= 0) {
         (self->ptrObj)->setConnectedConstraints(cellType, connectedLambda);
     }

     if (persistenceLambda >= 0 && persistenceDiffusion >= 0 && persistenceTime >= 0) {
         (self->ptrObj)->setPersistenceConstraints(cellType, persistenceLambda, 
                 persistenceTime, persistenceDiffusion);
     }

     if (otherCellType != -1 && adhesion != -1) {
         (self->ptrObj)->setAdhesionBetweenTypes(cellType, otherCellType, 
                 adhesion);
     }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_addCell(PyCpm2d* self, PyObject* args)
{
    int x;
    int y;
    int type;
    if (! PyArg_ParseTuple(args, "iii", &type, &x, &y))
        return Py_False;

    (self->ptrObj)->addCell(x, y, type);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_addCell(PyCpm3d* self, PyObject* args)
{
    int x;
    int y;
    int z;
    int type;
    if (! PyArg_ParseTuple(args, "iiii", &type, &x, &y, &z))
        return Py_False;

    (self->ptrObj)->addCell(x, y, z, type);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_setPoint(PyCpm2d* self, PyObject* args)
{
    int x;
    int y;
    int cellId;
    int type;
    if (! PyArg_ParseTuple(args, "iiii", &x, &y, &cellId, &type))
        return Py_False;

    (self->ptrObj)->setPoint(x, y, cellId, type);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_setPoint(PyCpm3d* self, PyObject* args)
{
    int x;
    int y;
    int z;
    int cellId;
    int type;
    if (! PyArg_ParseTuple(args, "iiii", &x, &y, &z, &cellId, &type))
        return Py_False;

    (self->ptrObj)->setPoint(x, y, z, cellId, type);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_updateType(PyCpm2d* self, PyObject* args)
{
    int id, type;

    if (! PyArg_ParseTuple(args, "ii", &id, &type))
        return Py_False;
    (self->ptrObj)->updateType(id, type);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_updateType(PyCpm2d* self, PyObject* args)
{
    int id, type;

    if (! PyArg_ParseTuple(args, "ii", &id, &type))
        return Py_False;
    (self->ptrObj)->updateType(id, type);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject * PyCpm2d_run(PyCpm2d* self, PyObject* args)
{
    int ticks;

    if (! PyArg_ParseTuple(args, "i", &ticks))
        return Py_False;
    (self->ptrObj)->run(ticks);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_run(PyCpm3d* self, PyObject* args)
{
    int ticks;

    if (! PyArg_ParseTuple(args, "i", &ticks))
        return Py_False;

    (self->ptrObj)->run(ticks);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_runAsync(PyCpm2d* self, PyObject* args)
{
    int ticks;

    if (! PyArg_ParseTuple(args, "i", &ticks))
        return Py_False;

    (self->ptrObj)->runAsync(ticks);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_runAsync(PyCpm3d* self, PyObject* args)
{
    int ticks;

    if (! PyArg_ParseTuple(args, "i", &ticks))
        return Py_False;

    (self->ptrObj)->runAsync(ticks);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_join(PyCpm3d* self, PyObject* args)
{
    (self->ptrObj)->join();

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_join(PyCpm3d* self, PyObject* args)
{
    (self->ptrObj)->join();

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_getField(PyCpm2d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {2,dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(3, shape, NPY_DOUBLE, 
            (self->ptrObj)->getField());
    return arr;
}


static PyObject * PyCpm3d_getField(PyCpm3d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {3,dimension, dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(4, shape, NPY_DOUBLE, 
            (self->ptrObj)->getField());
    return arr;
}

static PyObject * PyCpm2d_getState(PyCpm2d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(2, shape, NPY_INT, 
            (self->ptrObj)->getData());
    return arr;
}


static PyObject * PyCpm3d_getState(PyCpm3d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {dimension, dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(3, shape, NPY_INT, 
            (self->ptrObj)->getData());
    return arr;
}


static PyObject * PyCpm2d_getActState(PyCpm2d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(2, shape, NPY_INT, 
            (self->ptrObj)->getActData());
    return arr;
}

static PyObject * PyCpm3d_getActState(PyCpm3d* self, PyObject* args)
{
    int dimension = (self->ptrObj)->getDimension();
    npy_intp shape[] = {dimension, dimension, dimension};
    PyObject*  arr = PyArray_SimpleNewFromData(3, shape, NPY_INT, 
            (self->ptrObj)->getActData());
    return arr;
}


static PyObject * PyCpm2d_getCentroids(PyCpm2d* self, PyObject* args)
{
    auto centroids = (self->ptrObj)->getCentroids();

    npy_intp const dims[2] = {int(centroids.size()), 2};
    PyArrayObject* output = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
    if (!output) 
        return 0;
    double* data = (double*)output->data;
    for(int i = 0; i < centroids.size(); i++) {
        data[i*2 + 0] = centroids[i].x;
        data[i*2 + 1] = centroids[i].y;
    }
    return (PyObject*)output;
}

static PyObject * PyCpm3d_getCentroids(PyCpm3d* self, PyObject* args)
{
    auto centroids = (self->ptrObj)->getCentroids();

    npy_intp const dims[2] = {int(centroids.size()), 3};
    PyArrayObject* output = (PyArrayObject*) PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
    if (!output) 
        return 0;
    double* data = (double*)output->data;
    for(int i = 0; i < centroids.size(); i++) {
        data[i*3 + 0] = centroids[i].x;
        data[i*3 + 1] = centroids[i].y;
        data[i*3 + 2] = centroids[i].z;
    }
    return (PyObject*)output;
}

static PyObject * PyCpm2d_overwriteCell(PyCpm2d* self, PyObject* args)
{
    PyObject *arg=NULL;
    int type;
    if (!PyArg_ParseTuple(args, "Oi", &arg, &type)) return NULL;
    int dims = PyArray_NDIM(arg);
    npy_intp* dim_vals = PyArray_DIMS(arg);
    int nrOfPoints = dim_vals[0];
    int b = dim_vals[1];
    
    (self->ptrObj)->addCell(type);

    int dimension = (self->ptrObj)->getDimension();
    for (int i = 0; i < nrOfPoints; i++) {
       
        int x = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, 0));
        int y = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, 1));
        (self->ptrObj)->updatePoint(y, x, type);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_overwriteCell(PyCpm3d* self, PyObject* args)
{
    PyObject *arg=NULL;
    int type;
    if (!PyArg_ParseTuple(args, "Oi", &arg, &type)) return NULL;
    int dims = PyArray_NDIM(arg);
    npy_intp* dim_vals = PyArray_DIMS(arg);
    int nrOfPoints = dim_vals[0];
    //int y = dim_vals[1];
    
    (self->ptrObj)->addCell(type);

    int dimension = (self->ptrObj)->getDimension();
    for (int i = 0; i < nrOfPoints; i++) {
        int x = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, 0));
        int y = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, 1));
        int z = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, 2));
        (self->ptrObj)->updatePoint(z, y, x, type);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm2d_initializeFromArray(PyCpm2d* self, PyObject* args)
{
    PyObject *arg=NULL;
    int count;
    if (!PyArg_ParseTuple(args, "Oi", &arg, &count)) return NULL;
    int dims = PyArray_NDIM(arg);
    npy_intp* dim_vals = PyArray_DIMS(arg);
    int x = dim_vals[0];
    int y = dim_vals[1];

    int dimension = (self->ptrObj)->getDimension();

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            unsigned int cellId = *static_cast<unsigned int*>(PyArray_GETPTR2(arg, i, j));
            int id = cellId & 16777215U;
            int type = cellId >> 24;
            if (cellId != 0) {
                (self->ptrObj)->setPoint(j, i, id, type);
            }
        }
    }
    (self->ptrObj)->updateCellProps(count);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * PyCpm3d_initializeFromArray(PyCpm3d* self, PyObject* args)
{
    PyObject *arg=NULL;
    int count;
    if (!PyArg_ParseTuple(args, "Oi", &arg, &count)) return NULL;
    int dims = PyArray_NDIM(arg);
    npy_intp* dim_vals = PyArray_DIMS(arg);
    int x = dim_vals[0];
    int y = dim_vals[1];
    int z = dim_vals[2];

    int dimension = (self->ptrObj)->getDimension();


    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                long cellId = *static_cast<long*>(PyArray_GETPTR3(arg, i, j, k));
                int id = cellId & 16777215U;
                int type = cellId >> 24;
                if (cellId != 0) {
                    (self->ptrObj)->setPoint(k, j, i, id, type);
                }
            }
        }
    }
    (self->ptrObj)->updateCellProps(count);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PyCpm2d_methods[] = {
    { "initialize_from_array", (PyCFunction)PyCpm2d_initializeFromArray, METH_VARARGS, "initialize simulation from array" },
    { "add_cell", (PyCFunction)PyCpm2d_addCell, METH_VARARGS, "add new cell at location" },
    { "overwrite_cell", (PyCFunction)PyCpm2d_overwriteCell, METH_VARARGS, "add new cell at location" },
    { "set_point", (PyCFunction)PyCpm2d_setPoint, METH_VARARGS, "set point on lattice for cell that already exists" },
    { "run", (PyCFunction)PyCpm2d_run, METH_VARARGS, "run for certain number of ticks" },
    { "run_async", (PyCFunction)PyCpm2d_runAsync, METH_VARARGS, "run for certain number of ticks in seperate thread" },
    { "join", (PyCFunction)PyCpm2d_join, METH_VARARGS, "join if simulation is running asynchronously" },
    { "get_state", (PyCFunction)PyCpm2d_getState, METH_VARARGS, "get state of CPM lattice" },
    { "get_field", (PyCFunction)PyCpm2d_getField, METH_VARARGS, "get chemotaxis field of CPM" },
    { "get_act_state", (PyCFunction)PyCpm2d_getActState, METH_VARARGS, "get state of CPM act lattice" },
    { "get_centroids", (PyCFunction)PyCpm2d_getCentroids, METH_VARARGS, "get centroids of cells" },
    { "set_constraints", (PyCFunction)PyCpm2d_setConstraints, METH_VARARGS | METH_KEYWORDS, "get state of CPM lattice" },
    { "update_type", (PyCFunction)PyCpm2d_updateType, METH_VARARGS, "get state of CPM lattice" },
    {NULL}  /* Sentinel */
};

static PyMethodDef PyCpm3d_methods[] = {
    { "initialize_from_array", (PyCFunction)PyCpm3d_initializeFromArray, METH_VARARGS, "initialize simulation from array" },
    { "add_cell", (PyCFunction)PyCpm3d_addCell, METH_VARARGS, "add new cell at location" },
    { "overwrite_cell", (PyCFunction)PyCpm3d_overwriteCell, METH_VARARGS, "add new cell at location" },
    { "set_point", (PyCFunction)PyCpm3d_setPoint, METH_VARARGS, "set point on lattice for cell that already exists" },
    { "run", (PyCFunction)PyCpm3d_run, METH_VARARGS, "run for certain number of ticks" },
    { "run_async", (PyCFunction)PyCpm3d_runAsync, METH_VARARGS, "run for certain number of ticks in seperate thread" },
    { "join", (PyCFunction)PyCpm3d_join, METH_VARARGS, "join if simulation is running asynchronously" },
    { "get_state", (PyCFunction)PyCpm3d_getState, METH_VARARGS, "get state of CPM lattice" },
    { "get_field", (PyCFunction)PyCpm3d_getField, METH_VARARGS, "get chemotaxis field of CPM" },
    { "get_act_state", (PyCFunction)PyCpm3d_getActState, METH_VARARGS, "get state of CPM act lattice" },
    { "get_centroids", (PyCFunction)PyCpm3d_getCentroids, METH_VARARGS, "get centroids of cells" },
    { "set_constraints", (PyCFunction)PyCpm3d_setConstraints, METH_VARARGS | METH_KEYWORDS, "get state of CPM lattice" },
    { "update_type", (PyCFunction)PyCpm3d_updateType, METH_VARARGS, "get state of CPM lattice" },
    {NULL}  /* Sentinel */
};



PyMODINIT_FUNC PyInit_cpm(void)
{
    import_array();

    PyObject* m;

    PyCpm3dType.tp_new = PyType_GenericNew;
    PyCpm3dType.tp_basicsize=sizeof(PyCpm3d);
    PyCpm3dType.tp_dealloc=(destructor) PyCpm3d_dealloc;
    PyCpm3dType.tp_flags=Py_TPFLAGS_DEFAULT;
    PyCpm3dType.tp_doc="3d CPM objects";
    PyCpm3dType.tp_methods=PyCpm3d_methods;
    PyCpm3dType.tp_init=(initproc)PyCpm3d_init;

    if (PyType_Ready(&PyCpm3dType) < 0)
        return NULL;

    PyCpm2dType.tp_new = PyType_GenericNew;
    PyCpm2dType.tp_basicsize=sizeof(PyCpm2d);
    PyCpm2dType.tp_dealloc=(destructor) PyCpm2d_dealloc;
    PyCpm2dType.tp_flags=Py_TPFLAGS_DEFAULT;
    PyCpm2dType.tp_doc="2d CPM objects";
    PyCpm2dType.tp_methods=PyCpm2d_methods;
    PyCpm2dType.tp_init=(initproc)PyCpm2d_init;

    if (PyType_Ready(&PyCpm2dType) < 0)
        return NULL;
    
    m = PyModule_Create(&cpmmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&PyCpm3dType);
    PyModule_AddObject(m, "Cpm3d", (PyObject *)&PyCpm3dType); 
    Py_INCREF(&PyCpm2dType);
    PyModule_AddObject(m, "Cpm2d", (PyObject *)&PyCpm2dType); 
    return m;
}
