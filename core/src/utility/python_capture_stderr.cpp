//
// Copyright (C) 2011 Mateusz Loskot <mateusz@loskot.net>
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Blog article: http://mateusz.loskot.net/?p=2819

#include "python_capture_stderr.h"

#include <functional>
#include <iostream>
#include <string>
#include <Python.h>

namespace emb
{

typedef std::function<void(std::string)> stderr_write_type;

PyObject* Stderr_write(PyObject* self, PyObject* args)
{
    std::size_t written(0);
    Stderr* selfimpl = reinterpret_cast<Stderr*>(self);
    if (selfimpl->write)
    {
        char* data;
        if (!PyArg_ParseTuple(args, "s", &data))
            return 0;

        std::string str(data);
        selfimpl->write(str);
        written = str.size();
    }
    return PyLong_FromSize_t(written);
}

PyObject* Stderr_flush(PyObject* self, PyObject* args)
{
    // no-op
    return Py_BuildValue("");
}

PyMethodDef Stderr_methods[] =
{
    {"write", Stderr_write, METH_VARARGS, "sys.stderr.write"},
    {"flush", Stderr_flush, METH_VARARGS, "sys.stderr.write"},
    {0, 0, 0, 0} // sentinel
};

PyTypeObject StderrType =
{
    PyVarObject_HEAD_INIT(0, 0)
    "emb.StderrType",     /* tp_name */
    sizeof(Stderr),       /* tp_basicsize */
    0,                    /* tp_itemsize */
    0,                    /* tp_dealloc */
    0,                    /* tp_print */
    0,                    /* tp_getattr */
    0,                    /* tp_setattr */
    0,                    /* tp_reserved */
    0,                    /* tp_repr */
    0,                    /* tp_as_number */
    0,                    /* tp_as_sequence */
    0,                    /* tp_as_mapping */
    0,                    /* tp_hash  */
    0,                    /* tp_call */
    0,                    /* tp_str */
    0,                    /* tp_getattro */
    0,                    /* tp_setattro */
    0,                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,   /* tp_flags */
    "emb.Stderr objects", /* tp_doc */
    0,                    /* tp_traverse */
    0,                    /* tp_clear */
    0,                    /* tp_richcompare */
    0,                    /* tp_weaklistoffset */
    0,                    /* tp_iter */
    0,                    /* tp_iternext */
    Stderr_methods,       /* tp_methods */
    0,                    /* tp_members */
    0,                    /* tp_getset */
    0,                    /* tp_base */
    0,                    /* tp_dict */
    0,                    /* tp_descr_get */
    0,                    /* tp_descr_set */
    0,                    /* tp_dictoffset */
    0,                    /* tp_init */
    0,                    /* tp_alloc */
    0,                    /* tp_new */
};

PyModuleDef embmodule =
{
    PyModuleDef_HEAD_INIT,
    "emb", 0, -1, 0,
};

// Internal state
PyObject* g_stderr;
PyObject* g_stderr_saved;

PyMODINIT_FUNC PyInit_emb(void)
{
    g_stderr = 0;
    g_stderr_saved = 0;

    StderrType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&StderrType) < 0)
        return 0;

    PyObject* m = PyModule_Create(&embmodule);
    if (m)
    {
        Py_INCREF(&StderrType);
        PyModule_AddObject(m, "Stderr", reinterpret_cast<PyObject*>(&StderrType));
    }
    return m;
}

void set_stderr(stderr_write_type write)
{
    if (!g_stderr)
    {
        g_stderr_saved = PySys_GetObject("stderr"); // borrowed
        g_stderr = StderrType.tp_new(&StderrType, 0, 0);
    }

    Stderr* impl = reinterpret_cast<Stderr*>(g_stderr);
    impl->write = write;
    PySys_SetObject("stderr", g_stderr);
}

void reset_stderr()
{
    if (g_stderr_saved)
        PySys_SetObject("stderr", g_stderr_saved);

    Py_XDECREF(g_stderr);
    g_stderr = 0;
}

} // namespace emb

